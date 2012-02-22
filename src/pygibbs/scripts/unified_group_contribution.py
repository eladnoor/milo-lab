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

class UnknownReactionEnergyError(Exception):
    pass

class UnifiedGroupContribution(object):
    
    def __init__(self, db, html_writer):
        self.db = db
        self.html_writer = html_writer
    
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
        for i, col_indices in col_mapping.iteritems():
            self.b[0, i] = np.mean(b[0, col_indices])
            self.anchored[0, i] = anchored[0, col_indices].max()
            self.obs_ids.append(','.join([obs_collection.observations[j].obs_id for j in col_indices]))
            self.obs_types.append(','.join(set(obs_collection.observations[j].obs_type for j in col_indices)))
    
        n_groups = len(groups_data.GetGroupNames()) # number of groups
        self.G = np.matrix(np.zeros((len(self.cids), n_groups)))
        self.has_groupvec = np.matrix(np.zeros((len(self.cids), 1)))
        for i, cid in enumerate(self.cids):
            if cid2groupvec[cid] is not None:
                self.has_groupvec[i, 0] = 1
                self.G[i, :] = cid2groupvec[cid].Flatten()

    def ToDatabase(self):
        kegg = Kegg.getInstance()
        self.db.SaveSparseNumpyMatrix('ugc_S', self.S)
        self.db.SaveSparseNumpyMatrix('ugc_G', self.G)
        self.db.SaveNumpyMatrix('ugc_b', self.b.T)
        self.db.SaveNumpyMatrix('ugc_anchored', self.anchored.T)
        self.db.CreateTable('ugc_compounds', 'cid INT, name TEXT')
        for cid in self.cids:
            self.db.Insert('ugc_compounds', [cid, kegg.cid2name(cid)])
        self.db.CreateTable('ugc_observations', 'row INT, id TEXT, type TEXT')
        for i in xrange(len(self.obs_ids)):
            self.db.Insert('ugc_observations', [i, self.obs_ids[i], self.obs_types[i]])
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
        for rowdict in self.db.DictReader('ugc_observations'):
            self.obs_ids.append(rowdict['id'])
            self.obs_types.append(rowdict['type'])
    
    def NormalizeAnchors(self):
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
        
        # removed zero rows (compounds) from S
        used_cid_indices = list(np.nonzero(np.sum(abs(self.S), 1))[0].flat)
        self.S = self.S[used_cid_indices, :]
        self.G = self.G[used_cid_indices, :]
        self.cids = [self.cids[i] for i in used_cid_indices]
        self.has_groupvec = self.has_groupvec[used_cid_indices, :]

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
        est = np.matrix(np.zeros((3, r.shape[1]))) * np.nan
        
        try:
            # calculate the contributions of compounds
            g_S, PC_S, PL_S = UnifiedGroupContribution.regress(S, b)
            r_C = PC_S * r # the part of 'r' which is in the column-space of S
            r_L = PL_S * r # the part of 'r' which is in the left-null-space of S
            est[0, :] = g_S * r_C
        except np.linalg.linalg.LinAlgError:
            return est

        bad_compounds = np.where(self.has_groupvec == False)[0]
        reactions_with_groupvec = []
        for i in xrange(S.shape[1]):
            if np.all(abs(S[bad_compounds, i]) < 1e-10):
                reactions_with_groupvec.append(i)

        try:
            GS = self.G.T * S[:, reactions_with_groupvec]
            g_GS, PC_GS, PL_GS = UnifiedGroupContribution.regress(GS, b[:, reactions_with_groupvec])
            r_groupvec = self.G.T * r
            est[1, :] = g_GS * r_groupvec
            for k in xrange(r.shape[1]):
                if (abs(r[bad_compounds, k]) > 1e-10).any():
                    # this reaction involves compounds which don't have groupvectors
                    est[1, k] = np.nan
                elif (abs(PL_GS * r_groupvec[:, k]) > 1e-10).any():
                    # this reaction's groupvector is not in the column space of GS
                    est[1, k] = np.nan
        
        except np.linalg.linalg.LinAlgError:
            return est

        r_groupvec_L = self.G.T * r_L
        est[2, :] = g_GS * r_groupvec_L + est[0, :]
        for k in xrange(r.shape[1]):
            if np.any(abs(r_L[bad_compounds, k]) > 1e-10):
                # this reaction involves compounds which don't have groupvectors
                est[2, k] = np.nan
            elif np.any(abs(PL_GS * r_groupvec_L) > 1e-10):
                # this reaction's groupvector is not in the column space of GS
                est[2, k] = np.nan

            if np.any(abs(r_L[:, k]) > 1e-10):
                # this reaction is not in the column space of S 
                est[0, k] = np.nan
        
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
        sparse = dict((cids[c], S_row[c]) for c in active_cids)
        r = Reaction("", sparse)
        return r.FullReactionString(show_cids=False)

    def Report(self, est):
        legend = ['PRC', 'PGC', 'UGC']
        resid = est - np.repeat(self.b.T, est.shape[0], 1).T
        
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
            #rowdict['ID'] = self.obs_ids[i]
            #rowdict['type'] = self.obs_types[i]
            rowdict['reaction (anch)'] = UnifiedGroupContribution.row2hypertext(self.S_anchored[:, i], self.cids_anchored)
            rowdict['observed (anch)'] = "%.1f" % self.b_anchored[0, i]
            rowdict['reaction'] = UnifiedGroupContribution.row2hypertext(self.S[:, i], self.cids)
            rowdict['observed'] = "%.1f" % self.b[0, i]
            for j in xrange(est.shape[0]):
                rowdict[legend[j]] = "%.1f" % resid[j, i]
            rowdict['key'] = abs(resid[1, i])
            rowdicts.append(rowdict)
        rowdicts.sort(key=lambda x:x['key'], reverse=True)
        self.html_writer.write_table(rowdicts,
            #headers=['#', 'ID', 'type', 'reaction', 'observed'] + legend)
            headers=['#', 'reaction (anch)', 'observed (anch)', 'reaction', 'observed'] + legend)
        self.html_writer.div_end()
        
    def Loo(self):
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')

        n = self.S.shape[1]
        est = np.matrix(np.zeros((3, n))) * np.nan
        for i in xrange(n):
            if self.obs_types[i] == 'formation':
                continue
            if np.all(abs(self.S[:, i]) < 1e-10): # empty reaction
                continue
            no_i = range(0, i) + range(i+1, n)
            e = self.Estimate(self.S[:, no_i], self.b[:, no_i], self.S[:, i:i+1])
            est[:, i:i+1] = e
            print 'dG0\' = %7.1f' % self.b[0, i],
            print '| PRC = %5.1f' % abs(self.b[0, i] - est[0, i]),
            print '| PGC = %5.1f' % abs(self.b[0, i] - est[1, i]),
            print '| UGC = %5.1f' % abs(self.b[0, i] - est[2, i]),
            print '| %s' % UnifiedGroupContribution.row2string(self.S[:, i], self.cids)        
        self.Report(est)
        
    def Fit(self):
        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
        est = self.Estimate(self.S, self.b, self.S)
        self.Report(est)
    
if __name__ == "__main__":
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    html_writer = HtmlWriter('../res/ugc.html')
    ugc = UnifiedGroupContribution(db, html_writer)
    if False:
        ugc.LoadData()
        ugc.ToDatabase()
        sys.exit(0)
    else:
        ugc.FromDatabase()
    ugc.NormalizeAnchors()
    ugc.Fit()
    ugc.Loo()
    