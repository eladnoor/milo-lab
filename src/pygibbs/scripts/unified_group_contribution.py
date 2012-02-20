import numpy as np
from matplotlib.mlab import rms_flat
import matplotlib.pyplot as plt

from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.group_observation import GroupObervationCollection
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.kegg import Kegg

class UnknownReactionEnergyError(Exception):
    pass

class UnifiedGroupContribution(object):
    
    def __init__(self, db, html_writer):
        self.db = db
        self.html_writer = html_writer
        
    def LoadData(self, db):
        groups_data = GroupsData.FromDatabase(self.db, transformed=False)
    
        cid2nH_nMg = {}
        cid2error = {}
        cid2groupvec = {}
        
        for row in db.DictReader('pgc_groupvector'):
            cid = row['cid']
            gv_str = row['groupvec']
            if gv_str is not None:
                groupvec = GroupVector.FromJSONString(groups_data, gv_str)
            else:
                groupvec = None
            cid2groupvec[cid] = groupvec
            cid2error[cid] = row['err']
            cid2nH_nMg[cid] = (row['nH'], row['nMg'])
        
        obs_collection = GroupObervationCollection.FromDatabase(db=db,
                                table_name='pgc_observations', transformed=False)
        self.cids, self.S, self.b, self.anchored = obs_collection.GetStoichiometry()
        obs_ids = [obs.obs_id for obs in obs_collection.observations]
        obs_types = [obs.obs_type for obs in obs_collection.observations]
    
        n_groups = len(groups_data.GetGroupNames()) # number of groups
        G = np.zeros((len(self.cids), n_groups))
        self.has_groupvec = np.zeros((len(self.cids), 1))
        for i, cid in enumerate(self.cids):
            if cid2groupvec[cid] is not None:
                self.has_groupvec[i, 0] = 1
                G[i, :] = cid2groupvec[cid].Flatten()
        
        self.G = G
        self.obs_types = obs_types
        self.obs_ids = obs_ids
        
        self.S_anchored = np.zeros(self.S.shape)
        self.b_anchored = np.zeros(self.b.shape)
    
    @staticmethod
    def regress(A, y):
        U, s, V = np.linalg.svd(A, full_matrices=True)
        r = len(np.where(s > 1e-10)[0]) # the rank of A
        inv_S = np.zeros(A.shape).T
        for i in xrange(r):
            inv_S[i, i] = 1.0 / s[i]
        x = np.dot(np.dot(np.dot(y, V.T), inv_S), U.T)
        P_C = np.dot(U[:,:r], U[:,:r].T) # a projection matrix onto the row-space of A
        P_L = np.dot(U[:,r:], U[:,r:].T) # a projection matrix onto the null-space of A
        return x, P_C, P_L
    
    def Estimate(self, S, b, r):
        """
            Given
                - group matrix: G (m x g)
                - stoichiometric matrix: S (m x n)
                - Gibbs energies: b (1 x n) 
                - a new reaction: r (m x k)
            Return:
                The estimated Gibbs energy of 'r'
        """
        # calculate the contributions of compounds
        g_S, PC_S, PL_S = UnifiedGroupContribution.regress(S, b)
        r_C = np.dot(PC_S, r) # the part of 'r' which is in the column-space of S
        r_L = np.dot(PL_S, r) # the part of 'r' which is in the left-null-space of S
        
        est_prc = np.dot(g_S, r_C)
        for k in xrange(r.shape[1]):
            if np.any(abs(r_L[:, k]) > 1e-10):
                # this reaction is not in the column space of S 
                est_prc[0, k] = np.nan

        bad_compounds = np.where(self.has_groupvec==0)[1]
        reactions_to_use = []
        for i in xrange(S.shape[1]):
            if np.all(S[bad_compounds, i] == 0):
                reactions_to_use.append(i)
        GS = np.dot(self.G.T, S[:, reactions_to_use])
        g_GS, PC_GS, PL_GS = UnifiedGroupContribution.regress(GS, b)

        est_pgc = np.zeros((1, r.shape[1]))
        for k in xrange(r.shape[1]):
            if np.any(abs(r[bad_compounds, k]) > 1e-10):
                # this reaction involves compounds which don't have groupvectors
                est_pgc[0, k] = np.nan
                continue

            r_groupvec = np.dot(self.G.T, r[:, k])
            if np.any(abs(np.dot(PL_GS, r_groupvec)) > 1e-10):
                # this reaction's groupvector is not in the column space of GS
                est_pgc[0, k] = np.nan
                continue
            est_pgc[0, k] = np.dot(g_GS, r_groupvec)[0]

        est_ugc = np.zeros((1, r.shape[1]))
        for k in xrange(r.shape[1]):
            if np.any(abs(r_L[bad_compounds, k]) > 1e-10):
                # this reaction involves compounds which don't have groupvectors
                est_ugc[0, k] = np.nan
                continue

            r_groupvec_L = np.dot(self.G.T, r_L[:, k])
            if np.any(abs(np.dot(PL_GS, r_groupvec_L)) > 1e-10):
                # this reaction's groupvector is not in the column space of GS
                est_ugc[0, k] = np.nan
                continue

            est_ugc[0, k] = np.dot(g_GS, r_groupvec_L)[0] + est_prc[0, k]
            
        return np.vstack([est_prc, est_pgc, est_ugc])
    
    def Report(self, est):
        kegg = Kegg.getInstance()
        
        legend = ['PRC', 'PGC', 'UGC']
        resid = est - np.repeat(self.b, est.shape[0], 0)
        
        fig = plt.figure(figsize=(6,6))
        rms = []
        for j in xrange(est.shape[0]):
            used_reactions = np.where(np.isfinite(resid[j, :]))[0]
            rms.append(rms_flat(resid[j, used_reactions]))
            plt.plot(self.b[0, used_reactions],
                     resid[j, used_reactions],
                     '.', figure=fig, label=legend[j])
        
        self.html_writer.write_ul(["RMSE(%s) = %.1f" % (legend[j], rms[j])
                                   for j in xrange(est.shape[0])])
        
        plt.xlabel(r"$\Delta_r G^{'\circ}$ observed [kJ/mol]")
        plt.ylabel(r"$\Delta_r G^{'\circ}$ residual [kJ/mol]")
        plt.legend()
        self.html_writer.embed_matplotlib_figure(fig)
        
        self.html_writer.insert_toggle(start_here=True, label="Show table")
        rowdicts = []
        for i in used_reactions:
            rowdict = {}
            S_row = self.S[:, i] + self.S_anchored[:, i]
            sparse = dict((self.cids[c], S_row[c]) for c in np.nonzero(S_row)[0])
            rowdict['reaction'] = kegg.sparse_to_hypertext(sparse, show_cids=False)
            rowdict['ID'] = self.obs_ids[i]
            rowdict['type'] = self.obs_types[i]
            rowdict['observed'] = "%.1f" % (self.b[0, i] + self.b_anchored[0, i])
            for j in xrange(est.shape[0]):
                rowdict[legend[j]] = "%.1f" % resid[j, i]
            rowdict['key'] = abs(resid[1, i])
            rowdicts.append(rowdict)
        rowdicts.sort(key=lambda x:x['key'], reverse=True)
        self.html_writer.write_table(rowdicts,
            headers=['#', 'ID', 'type', 'reaction', 'observed'] + legend)
        self.html_writer.div_end()
        
    def Loo(self):
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')

        n = self.S.shape[1]
        est = np.zeros((3, n))
        for i in xrange(n):
            no_i = range(0, i) + range(i+1, n)
            e = self.Estimate(self.S[:, no_i], self.b[:, no_i], self.S[:, i:i+1])
            est[:, i:i+1] = e
            print 'obs = %.1f' % self.b[0, i]
            print 'PGC =', self.b[0, i] - est[0, i]
            print 'PRC =', self.b[0, i] - est[1, i]
            print 'UGC =', self.b[0, i] - est[2, i]
        
        self.Report(est)
        
    def Fit(self):
        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
        est = self.Estimate(self.S, self.b, self.S)
        self.Report(est)
    
    def NormalizeAnchors(self):
        # now remove anchored data from S and leave only the data which will be 
        # used for calculating the group contributions
        anchored_cols = np.where(self.anchored==1)[1]
        g, P_C, P_L = UnifiedGroupContribution.regress(self.S[:, anchored_cols],
                                                       self.b[:, anchored_cols])

        # calculate the matrix and observations which are explained
        # by the anchored reactions
        self.S_anchored = np.dot(P_C, self.S)
        self.b_anchored = np.dot(g, self.S_anchored)
        
        # calculate the matrix and observations which are in the null-space
        # of the anchored reactions. in other words, b_L are the residuals
        # of the linear regression 
        self.S -= self.S_anchored
        self.b -= self.b_anchored
        
        # set epsilon-small values to absolute 0
        self.S[np.where(abs(self.S) < 1e-10)] = 0
        
if __name__ == "__main__":
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    html_writer = HtmlWriter('../res/ugc.html')
    ugc = UnifiedGroupContribution(db, html_writer)
    ugc.LoadData(db)
    ugc.NormalizeAnchors()
    ugc.Fit()
    ugc.Loo()
    