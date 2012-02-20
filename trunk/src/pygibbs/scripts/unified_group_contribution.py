import numpy as np
from matplotlib.mlab import rms_flat
import matplotlib.pyplot as plt

from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.group_observation import GroupObervationCollection
from toolbox.linear_regression import LinearRegression
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
        bad_compounds = np.where(self.has_groupvec==0)[1]
        
        est_pgc = np.zeros((1, r.shape[1])) * np.nan
        est_ugc = np.zeros((1, r.shape[1])) * np.nan

        # calculate the contributions of compounds
        g_S, PC_S, PL_S = UnifiedGroupContribution.regress(S, b)
        
        for k in xrange(r.shape[1]):
            if np.all(r[:, k] == 0):
                est_pgc[0, k] = 0
                est_ugc[0, k] = 0
                continue

            if np.any(r[bad_compounds, k] != 0):
                # cannot use group contribution for these reactions. Use only PRC
                
                if np.any(abs(np.dot(PL_S, r[:, k])) > 1e-10):
                    continue
                
                est_ugc[0, k] = np.dot(g_S, r[:, k])[0]
                
                
            if np.any(abs(np.dot(PL_S, r[:, k])) > 1e-10):
                # this reaction is not in the space of S, but we might still
                # be able to use group contribution
                continue
                
        bad_cols = np.where(self.has_groupvec==0)[1]
        if np.any(r[:, bad_cols] != 0):
            # cannot use group contribution for these reactions. Use only PRC
            
            # calculate the contributions of compounds
            g_S, P_C, P_L = UnifiedGroupContribution.regress(S, b)
            
            for i in np.where(np.sum(abs(np.dot(P_L, r)), 0) > 1e-10):
                est_ugc[0, i] = np.nan
        else:
            cont_GS, P_C1, P_L1 = UnifiedGroupContribution.regress(np.dot(self.G.T, S), b)
        
            # we must check that r is in the column-space of GS,
            # otherwise there is no way to estimate anything
            r_in_groups = np.dot(self.G.T, r)
            if np.any(abs(np.dot(P_L1, r_in_groups)) > 1e-10):
                raise UnknownReactionEnergyError("The given reaction is not in the group space of observed reaction")
            
            # convert the group contributions to compound contributions
            cont_G = np.dot(cont_GS, self.G.T)
            est_pgc = np.dot(cont_G, r)   # the full Gibbs energy attributed to Groups
            
            # calculate the contributions of compounds back-calculated from the gruops
            cont_S, P_C2, P_L2 = UnifiedGroupContribution.regress(S, b)
            
            r_C = np.dot(P_C2, r) # the part of 'r' which is in the column-space of S
            r_L = np.dot(P_L2, r) # the part of 'r' which is in the left-null-space of S
            
            g_C = np.dot(cont_S, r_C) # the Gibbs energy attributed to NIST directly
            g_L = np.dot(cont_G, r_L) # the residual Gibbs energy attributed to Groups
            est_ugc = g_C + g_L
        
        return est_pgc, est_ugc
    
    def Report(self, b_pgc, b_ugc):
        kegg = Kegg.getInstance()
        used_reactions = np.where(np.isfinite(b_pgc))[1]
        
        r_pgc = rms_flat(self.b[0, used_reactions] - b_pgc[0, used_reactions])
        r_ugc = rms_flat(self.b[0, used_reactions] - b_ugc[0, used_reactions])
        self.html_writer.write_ul(["RMSE(PGC) = %.1f" % r_pgc,
                                   "RMSE(UGC) = %.1f" % r_ugc])
        
        fig = plt.figure(figsize=(6,6))
        resid = np.vstack([self.b - b_pgc, self.b - b_ugc])
        plt.plot(self.b[:, used_reactions].T, resid[:, used_reactions].T, '.', figure=fig)
        plt.xlabel(r"$\Delta_r G^{'\circ}$ observed [kJ/mol]")
        plt.ylabel(r"$\Delta_r G^{'\circ}$ residual [kJ/mol]")
        plt.legend(['PGC', 'UGC'])
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
            rowdict['r<sub>PGC</sub>'] = "%.1f" % resid[0, i]
            rowdict['r<sub>UGC</sub>'] = "%.1f" % resid[1, i]
            rowdict['key'] = abs(resid[1, i])
            rowdicts.append(rowdict)
        rowdicts.sort(key=lambda x:x['key'], reverse=True)
        self.html_writer.write_table(rowdicts,
            headers=['#', 'ID', 'type', 'reaction', 'observed',
                     'r<sub>PGC</sub>', 'r<sub>UGC</sub>'])
        self.html_writer.div_end()
        
    def Loo(self):
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')

        n = self.S.shape[1]
        b_pgc = np.zeros(self.b.shape)
        b_ugc = np.zeros(self.b.shape)
        for i in xrange(n):
            no_i = range(0, i) + range(i+1, n)
            try:
                est_pgc, est_ugc = self..Estimate(self.S[:, no_i], self.b[:, no_i], self.S[:, i])
                b_pgc[0, i] = est_pgc[0]
                b_ugc[0, i] = est_ugc[0]
                print "%6.1f   %6.1f   %6.1f" % (self.b[0, i],
                                                 abs(self.b[0, i] - b_pgc[0, i]),
                                                 abs(self.b[0, i] - b_ugc[0, i]))
            except UnknownReactionEnergyError as e:
                print str(e)
                b_pgc[0, i] = np.nan
                b_ugc[0, i] = np.nan
        
        self.Report(b_pgc, b_ugc)
        
    def Fit(self):
        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
        b_pgc, b_ugc = self.Estimate(self.S, self.b, self.S)
        self.Report(b_pgc, b_ugc)
    
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
    