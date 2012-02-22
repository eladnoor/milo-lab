import numpy as np
from matplotlib.mlab import rms_flat
import matplotlib.pyplot as plt

from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.group_observation import GroupObervationCollection
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.kegg import Kegg
from toolbox.linear_regression import LinearRegression
import sys
from pygibbs.kegg_reaction import Reaction
from pygibbs.dissociation_constants import DissociationConstants

class UnknownReactionEnergyError(Exception):
    pass

class UnifiedGroupContribution(object):
    
    def __init__(self, db, html_writer):
        self.db = db
        self.html_writer = html_writer
        self.kegg = Kegg.getInstance()
    
    @staticmethod
    def regress(A, y):
        """
            Solves for x in the minimization of ||xA - y||
            using linear regression.
            
            Returns:
                x   - the regression result
                P_C - a projection matrix onto the column-space of A
                P_L - a projection matrix onto the left-null-space of A
                
            Note:
                Using x * r for extrapolating the values of 'y' to a new 
                column will only be valid if r is in the column-space of A.
                To check that, one must see that P_L * r == 0.
        """
        U, s, V = np.linalg.svd(A, full_matrices=True)
        r = len(np.where(s > 1e-10)[0]) # the rank of A
        inv_S = np.matrix(np.zeros(A.shape)).T
        for i in xrange(r):
            inv_S[i, i] = 1.0 / s[i]
        x = y * V.T * inv_S * U.T
        P_C = U[:,:r] * U[:,:r].T # a projection matrix onto the row-space of A
        P_L = U[:,r:] * U[:,r:].T # a projection matrix onto the null-space of A
        return x, P_C, P_L
    
    def LoadData(self):
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
        
        obs_collection = GroupObervationCollection.FromDatabase(db=self.db,
                                table_name='pgc_observations', transformed=False)
        self.cids, S, b, anchored = obs_collection.GetStoichiometry()
        self.S, col_mapping = LinearRegression.ColumnUnique(S)
        self.b = np.matrix(np.zeros((1, len(col_mapping)), dtype='float'))
        self.anchored = np.matrix(np.zeros((1, len(col_mapping)), dtype='int'))
        self.obs_ids = []
        self.obs_types = []
        self.obs_urls = []
        for i, col_indices in col_mapping.iteritems():
            self.b[0, i] = np.mean(b[0, col_indices])
            self.anchored[0, i] = anchored[0, col_indices].max()
            obs_list = [obs_collection.observations[j] for j in col_indices]
            self.obs_ids.append(', '.join([obs.obs_id for obs in obs_list]))
            self.obs_types.append(', '.join([obs.obs_type for obs in obs_list]))
            self.obs_urls.append(', '.join([obs.url for obs in obs_list]))
    
        n_groups = len(groups_data.GetGroupNames()) # number of groups
        self.G = np.matrix(np.zeros((len(self.cids), n_groups)))
        self.has_groupvec = np.matrix(np.zeros((len(self.cids), 1)))
        for i, cid in enumerate(self.cids):
            if cid2groupvec[cid] is not None:
                self.has_groupvec[i, 0] = 1
                self.G[i, :] = cid2groupvec[cid].Flatten()

    def ToDatabase(self):
        self.db.SaveSparseNumpyMatrix('ugc_S', self.S)
        self.db.SaveSparseNumpyMatrix('ugc_G', self.G)
        self.db.SaveNumpyMatrix('ugc_b', self.b.T)
        self.db.SaveNumpyMatrix('ugc_anchored', self.anchored.T)
        self.db.CreateTable('ugc_compounds', 'cid INT, name TEXT')
        for cid in self.cids:
            self.db.Insert('ugc_compounds', [cid, self.kegg.cid2name(cid)])
        self.db.CreateTable('ugc_observations', 'row INT, id TEXT, type TEXT, url TEXT')
        for i in xrange(len(self.obs_ids)):
            self.db.Insert('ugc_observations', [i, self.obs_ids[i], self.obs_types[i], self.obs_urls[i]])
        self.db.Commit()

    def FromDatabase(self):
        self.S = self.db.LoadSparseNumpyMatrix('ugc_S')
        self.G = self.db.LoadSparseNumpyMatrix('ugc_G')
        self.b = self.db.LoadNumpyMatrix('ugc_b').T
        self.anchored = self.db.LoadNumpyMatrix('ugc_anchored').T
        self.has_groupvec = np.sum(self.G, 1) > 0
        self.cids = []
        for rowdict in self.db.DictReader('ugc_compounds'):
            self.cids.append(int(rowdict['cid']))
        self.obs_ids = []
        self.obs_types = []
        self.obs_urls = []
        for rowdict in self.db.DictReader('ugc_observations'):
            self.obs_ids.append(rowdict['id'])
            self.obs_types.append(rowdict['type'])
            self.obs_urls.append(rowdict['url'])
    
    def SqueezeData(self, normalize_anchors=True):
        if normalize_anchors:
            # now remove anchored data from S and leave only the data which will be 
            # used for calculating the group contributions
            anchored_cols = list(np.where(self.anchored==1)[1].flat)
    
            g, P_C, P_L = UnifiedGroupContribution.regress(self.S[:, anchored_cols],
                                                           self.b[:, anchored_cols])
    
            # calculate the matrix and observations which are explained
            # by the anchored reactions
            self.S_anchored = P_C * self.S
            self.b_anchored = g * self.S_anchored
            self.cids_anchored = list(self.cids)
            
            # calculate the matrix and observations which are in the null-space
            # of the anchored reactions. in other words, b_L are the residuals
            # of the linear regression 
            self.S -= self.S_anchored
            self.b -= self.b_anchored
            
            # set epsilon-small values to absolute 0
            self.S[np.where(abs(self.S) < 1e-10)] = 0
        else:
            self.S_anchored = np.matrix(np.zeros(self.S.shape))
            self.b_anchored = np.matrix(np.zeros(self.b.shape))
            self.cids_anchored = list(self.cids)
        
        # removed zero rows (compounds) from S
        used_cid_indices = list(np.nonzero(np.sum(abs(self.S), 1))[0].flat)
        self.S = self.S[used_cid_indices, :]
        self.G = self.G[used_cid_indices, :]
        self.cids = [self.cids[i] for i in used_cid_indices]
        self.has_groupvec = self.has_groupvec[used_cid_indices, :]

    def Estimate(self, S, b, r, verbose=False):
        """
            Given
                - group matrix: G (m x g)
                - stoichiometric matrix: S (m x n)
                - Gibbs energies: b (1 x n) 
                - a new reaction: r (m x k)
            Return:
                The estimated Gibbs energy of 'r'
        """
        est = np.matrix(np.zeros((3, r.shape[1]))) * np.nan
        
        bad_compounds = np.where(self.has_groupvec == False)[0]
        reactions_with_groupvec = []
        for i in xrange(S.shape[1]):
            if np.all(abs(S[bad_compounds, i]) < 1e-10):
                reactions_with_groupvec.append(i)
        GS = self.G.T * S[:, reactions_with_groupvec]

        try:
            # calculate the contributions of compounds
            g_S, PC_S, PL_S = UnifiedGroupContribution.regress(S, b)
            g_GS, PC_GS, PL_GS = UnifiedGroupContribution.regress(GS, b[:, reactions_with_groupvec])
        except np.linalg.linalg.LinAlgError:
            return est

        r_C = PC_S * r # the part of 'r' which is in the column-space of S
        r_L = PL_S * r # the part of 'r' which is in the left-null-space of S


        r_groupvec = self.G.T * r     # a group-vector representation of r
        r_groupvec_L = self.G.T * r_L # a group-vector representation of r_L

        est[0, :] = g_S * r_C
        est[1, :] = g_GS * r_groupvec
        est[2, :] = g_GS * r_groupvec_L + est[0, :]
        
        for k in xrange(r.shape[1]):
            if verbose:
                print "r =", UnifiedGroupContribution.row2string(r[:, k], self.cids)
                print "r_C =", UnifiedGroupContribution.row2string(r_C[:, k], self.cids)
                print "r_L =", UnifiedGroupContribution.row2string(r_L[:, k], self.cids)

            if np.any(abs(r_L[:, k]) > 1e-10):
                # this reaction is not in the column space of S 
                est[0, k] = np.nan

            if (abs(r[bad_compounds, k]) > 1e-10).any():
                # this reaction involves compounds which don't have groupvectors
                est[1, k] = np.nan
            elif (abs(PL_GS * r_groupvec[:, k]) > 1e-10).any():
                # this reaction's groupvector is not in the column space of GS
                est[1, k] = np.nan

            if (abs(r_L[bad_compounds, k]) > 1e-10).any():
                # this reaction involves compounds which don't have groupvectors
                est[2, k] = np.nan
            elif (abs(PL_GS * r_groupvec_L[:, k]) > 1e-10).any():
                # this reaction's groupvector is not in the column space of GS
                est[2, k] = np.nan
        
        return est
    
    @staticmethod
    def row2hypertext(S_row, cids):
        kegg = Kegg.getInstance()
        active_cids = list(np.nonzero(S_row)[0].flat)
        sparse = dict((cids[c], S_row[c]) for c in active_cids)
        return kegg.sparse_to_hypertext(sparse, show_cids=False)

    @staticmethod
    def row2string(S_row, cids):
        active_cids = list(np.nonzero(S_row)[0].flat)
        sparse = dict((cids[c], S_row[c]) 
                      for c in active_cids 
                      if abs(S_row[c]) > 1e-10)
        r = Reaction("", sparse)
        return r.FullReactionString(show_cids=False)

    def Report(self, resid):
        legend = ['PRC', 'PGC', 'UGC']
        
        fig = plt.figure(figsize=(6,6))
        rms_list = []
        for j, label in enumerate(legend):
            used_reactions = list(np.where(np.isfinite(resid[j, :]))[1].flat)
            x = list(self.b[0, used_reactions].flat)
            y = list(resid[j, used_reactions].flat)
            rms = rms_flat(y)
            rms_list.append("RMSE(%s) = %.1f, N = %d" % (label, rms, len(x)))
            plt.plot(x, y, '.', figure=fig, label=label)
        
        self.html_writer.write_ul(rms_list)
        
        plt.xlabel(r"$\Delta_r G^{'\circ}$ observed [kJ/mol]")
        plt.ylabel(r"$\Delta_r G^{'\circ}$ residual [kJ/mol]")
        plt.legend()
        self.html_writer.embed_matplotlib_figure(fig)
        
        self.html_writer.insert_toggle(start_here=True, label="Show table")
        rowdicts = []
        for i in used_reactions:
            rowdict = {}
            rowdict['ID'] = self.obs_ids[i]
            rowdict['type'] = self.obs_types[i]
            rowdict['link'] = self.obs_urls[i]
            rowdict['reaction (anch)'] = UnifiedGroupContribution.row2hypertext(self.S_anchored[:, i], self.cids_anchored)
            rowdict['observed (anch)'] = "%.1f" % self.b_anchored[0, i]
            rowdict['reaction'] = UnifiedGroupContribution.row2hypertext(self.S[:, i], self.cids)
            rowdict['observed'] = "%.1f" % self.b[0, i]
            for j, label in enumerate(legend):
                rowdict[label] = "%.1f" % resid[j, i]
            rowdict['key'] = abs(resid[1, i])
            rowdicts.append(rowdict)
        #rowdicts.sort(key=lambda x:x['key'], reverse=True)
        self.html_writer.write_table(rowdicts,
            headers=['#', 'ID', 'type', 'reaction', 'observed'] + legend)
            #headers=['#', 'reaction (anch)', 'observed (anch)', 'reaction', 'observed'] + legend)
        self.html_writer.div_end()
        
    def Loo(self):
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')

        n = self.S.shape[1]
        resid = np.matrix(np.zeros((3, n))) * np.nan
        for i in xrange(n):
            if self.obs_types[i] == 'formation':
                continue
            if np.all(abs(self.S[:, i]) < 1e-10): # empty reaction
                continue
            no_i = range(0, i) + range(i+1, n)
            est = self.Estimate(self.S[:, no_i], self.b[:, no_i], 
                                self.S[:, i:i+1], verbose=True)
            resid[:, i:i+1] = est - self.b[0, i]
            print 'dG0\' = %7.1f' % self.b[0, i],
            print '| PRC = %5.1f' % resid[0, i],
            print '| PGC = %5.1f' % resid[1, i],
            print '| UGC = %5.1f' % resid[2, i],
            print '| %s' % UnifiedGroupContribution.row2string(self.S[:, i], self.cids)        
        self.Report(resid)
        
    def Fit(self):
        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
        est = self.Estimate(self.S, self.b, self.S)
        resid = est - np.repeat(self.b.T, est.shape[0], 1).T
        self.Report(resid)
    
    def Temp(self):
        n = self.S.shape[1]
        i = 81
        no_i = range(0, i) + range(i+1, n)
        S = self.S[:, no_i]
        b = self.b[:, no_i]
        r = self.S[:, i:i+1]
        g_S, PC_S, PL_S = UnifiedGroupContribution.regress(S, b)
        r_C = PC_S * r
        r_L = PL_S * r
        print "r =", UnifiedGroupContribution.row2string(r, self.cids)
        print "r_C =", UnifiedGroupContribution.row2string(r_C, self.cids)
        print "r_L =", UnifiedGroupContribution.row2string(r_L, self.cids)
    
if __name__ == "__main__":
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    html_writer = HtmlWriter('../res/ugc.html')
    
    if False: # reread observations from files
        dissociation = DissociationConstants.FromPublicDB()
        obs_collection = GroupObervationCollection.FromFiles(
                            html_writer=html_writer, 
                            dissociation=dissociation,
                            transformed=False)
        obs_collection.ToDatabase(db, 'pgc_observations')
        sys.exit(0)

    ugc = UnifiedGroupContribution(db, html_writer)
    if False:
        ugc.LoadData()
        ugc.ToDatabase()
        sys.exit(0)
    else:
        ugc.FromDatabase()
    ugc.SqueezeData(normalize_anchors=True)
    ugc.Fit()
    ugc.Loo()
    