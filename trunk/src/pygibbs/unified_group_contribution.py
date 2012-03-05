import numpy as np
from matplotlib.mlab import rms_flat
import matplotlib.pyplot as plt

from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.group_observation import GroupObervationCollection
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from pygibbs.kegg import Kegg
from toolbox.linear_regression import LinearRegression
import sys
from pygibbs.kegg_reaction import Reaction
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from optparse import OptionParser
from toolbox import util
import logging
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from toolbox.plotting import cdf
from pygibbs.pseudoisomer import PseudoisomerMap
import json
from pygibbs.thermodynamic_constants import default_T, R

class UnknownReactionEnergyError(Exception):
    pass

class UnifiedGroupContribution(PsuedoisomerTableThermodynamics):
    
    def __init__(self, db, html_writer=None, dissociation=None):
        PsuedoisomerTableThermodynamics.__init__(self, name="Unified Group Contribution")
        self.db = db
        self.html_writer = html_writer or NullHtmlWriter()
        self.dissociation = dissociation
        self.transformed = False
        self.epsilon = 1e-10
        self.kegg = Kegg.getInstance()
        
        self.STOICHIOMETRIC_TABLE_NAME = 'ugc_S'
        self.GROUP_TABLE_NAME = 'ugc_G'
        self.GIBBS_ENERGY_TABLE_NAME = 'ugc_b'
        self.ANCHORED_TABLE_NAME = 'ugc_anchored'
        self.COMPOUND_TABLE_NAME = 'ugc_compounds'
        self.OBSERVATION_TABLE_NAME = 'ugc_observations'
        self.GROUPVEC_TABLE_NAME = 'ugc_groupvectors'
        self.UNIQUE_OBSERVATION_TABLE_NAME = 'ugc_unique_observations'
        self.THERMODYNAMICS_TABLE_NAME = 'ugc_pseudoisomers'
        self.ERRORS_TABLE_NAME = 'ugc_errors'
        self.CONSERVATIONS_TABLE_NAME = 'ugc_conservations'

    def GetDissociationConstants(self):
        """
            Since loading the pKas takes time, this function is a lazy initialization
            of self.dissociation.
        """
        if self.dissociation is None:
            self.dissociation = DissociationConstants.FromPublicDB()
        return self.dissociation

    def LoadGroups(self, FromDatabase=False):
        #if self.transformed:
        #    fname = "../data/thermodynamics/groups_species_transformed.csv"
        #else:
        if FromDatabase and self.db.DoesTableExist('groups'):
            self.groups_data = GroupsData.FromDatabase(self.db,
                                                       transformed=self.transformed)
            self.group_decomposer = GroupDecomposer(self.groups_data)
        else:
            fname = "../data/thermodynamics/groups_species.csv"
            self.groups_data = GroupsData.FromGroupsFile(fname,
                                                         transformed=self.transformed)
            self.groups_data.ToDatabase(self.db)
            self.group_decomposer = GroupDecomposer(self.groups_data)
    
    def LoadObservations(self, FromDatabase=False):
        if FromDatabase and self.db.DoesTableExist(self.OBSERVATION_TABLE_NAME):
            logging.info("Reading observations from database")
            self.obs_collection = GroupObervationCollection.FromDatabase(
                                    db=self.db,
                                    table_name=self.OBSERVATION_TABLE_NAME,
                                    transformed=self.transformed)
        else:
            logging.info("Reading observations from files")
            dissociation = self.GetDissociationConstants()
            self.obs_collection = GroupObervationCollection.FromFiles(
                                    html_writer=self.html_writer, 
                                    dissociation=dissociation,
                                    transformed=self.transformed)
            self.obs_collection.ToDatabase(self.db, self.OBSERVATION_TABLE_NAME)
        
        self.obs_collection.ReportToHTML()

    def LoadGroupVectors(self, FromDatabase=False):
        self.cid2groupvec = {}
        self.cid2error = {}            

        if FromDatabase and self.db.DoesTableExist(self.GROUPVEC_TABLE_NAME):
            logging.info("Reading group-vectors from database")
            self.cid2nH_nMg = {}
            for row in self.db.DictReader(self.GROUPVEC_TABLE_NAME):
                cid = row['cid']
                gv_str = row['groupvec']
                if gv_str is not None:
                    groupvec = GroupVector.FromJSONString(self.groups_data,
                                                          gv_str)
                else:
                    groupvec = None
                self.cid2groupvec[cid] = groupvec
                self.cid2error[cid] = row['err']
                self.cid2nH_nMg[cid] = (row['nH'], row['nMg'])
        else:
            logging.info("Decomposing all compounds and calculating group vectors")
            self.html_writer.write('</br><b>All Groupvectors</b>\n')
            self.html_writer.insert_toggle(start_here=True)
            
            dissociation = self.GetDissociationConstants()

            # The group vector's pseudoisomers must be consistent with the
            # psuedoisomers used for the reverse transform.
            # Here we simply use the dictionary self.cid2nH_nMg from GroupObervationCollection
            self.cid2nH_nMg = self.obs_collection.cid2nH_nMg

            for cid in sorted(self.kegg.get_all_cids()):
                self.cid2groupvec[cid] = None
                self.cid2error[cid] = None
                if cid not in self.cid2nH_nMg:
                    self.cid2error[cid] = "Does not have data about major pseudoisomer"
                    continue
                nH, nMg = self.cid2nH_nMg[cid]
                diss_table = dissociation.GetDissociationTable(cid, False)
                if diss_table is None:
                    self.cid2error[cid] = "Does not have pKa data"
                    continue
                mol = diss_table.GetMol(nH=nH, nMg=nMg)
                if mol is None:
                    self.cid2error[cid] = "Does not have structural data"
                    continue
                try:
                    mol.RemoveHydrogens()
                    decomposition = self.group_decomposer.Decompose(mol, 
                                        ignore_protonations=False, strict=True)
                except GroupDecompositionError:
                    self.cid2error[cid] = "Could not be decomposed"
                    continue
                groupvec = decomposition.AsVector()
                if nH != groupvec.Hydrogens() or nMg != groupvec.Magnesiums():
                    err_msg = "C%05d's most abundant pseudoisomer is [nH=%d, nMg=%d], " \
                        "but the decomposition has [nH=%d, nMg=%d]. Skipping..." \
                        "" % (cid, nH, nMg, groupvec.Hydrogens(), groupvec.Magnesiums())
                    self.html_writer.write('</br>ERROR: %s\n' % err_msg)
                    self.cid2error[cid] = err_msg
                else:
                    self.cid2groupvec[cid] = groupvec
    
            self.db.CreateTable(self.GROUPVEC_TABLE_NAME,
                "cid INT, nH INT, nMg INT, groupvec TEXT, err TEXT")
            for cid, gv in sorted(self.cid2groupvec.iteritems()):
                nH, nMg = self.cid2nH_nMg.get(cid, (0, 0))
                if gv is not None:
                    gv_str = gv.ToJSONString()
                else:
                    gv_str = None
                err = self.cid2error[cid]
                self.db.Insert(self.GROUPVEC_TABLE_NAME,
                               [cid, nH, nMg, gv_str, err])
            self.db.Commit()
            self.html_writer.div_end()

    def _GenerateGroupMatrix(self, cids):
        n_groups = len(self.groups_data.GetGroupNames()) # number of groups
        G = np.matrix(np.zeros((len(cids), n_groups)))
        has_groupvec = np.matrix(np.zeros((len(cids), 1)))
        for i, cid in enumerate(cids):
            if self.cid2groupvec[cid] is not None:
                has_groupvec[i, 0] = 1
                G[i, :] = self.cid2groupvec[cid].Flatten()
        return G, has_groupvec

    def LoadData(self, FromDatabase=False):
        if FromDatabase and self.db.DoesTableExist(self.STOICHIOMETRIC_TABLE_NAME):
            self.S = self.db.LoadSparseNumpyMatrix(self.STOICHIOMETRIC_TABLE_NAME)
            self.G = self.db.LoadSparseNumpyMatrix(self.GROUP_TABLE_NAME)
            self.b = self.db.LoadNumpyMatrix(self.GIBBS_ENERGY_TABLE_NAME).T
            self.anchored = self.db.LoadNumpyMatrix(self.ANCHORED_TABLE_NAME).T
            self.has_groupvec = np.sum(self.G, 1) > 0
            self.cids = []
            for rowdict in self.db.DictReader(self.COMPOUND_TABLE_NAME):
                self.cids.append(int(rowdict['cid']))
            self.obs_ids = []
            self.obs_types = []
            self.obs_urls = []
            for rowdict in self.db.DictReader(self.OBSERVATION_TABLE_NAME):
                self.obs_ids.append(rowdict['id'])
                self.obs_types.append(rowdict['type'])
                self.obs_urls.append(rowdict['url'])
        else:
            self.cids, S, b, anchored = self.obs_collection.GetStoichiometry()
            self.S, col_mapping = LinearRegression.ColumnUnique(S)
            self.b = np.matrix(np.zeros((1, len(col_mapping)), dtype='float'))
            self.anchored = np.matrix(np.zeros((1, len(col_mapping)), dtype='int'))
            self.obs_ids = []
            self.obs_types = []
            self.obs_urls = []
            for i, col_indices in col_mapping.iteritems():
                self.b[0, i] = np.mean(b[0, col_indices])
                self.anchored[0, i] = anchored[0, col_indices].max()
                obs_list = [self.obs_collection.observations[j] for j in col_indices]
                self.obs_ids.append(', '.join([obs.obs_id for obs in obs_list]))
                self.obs_types.append(', '.join([obs.obs_type for obs in obs_list]))
                self.obs_urls.append(', '.join([obs.url for obs in obs_list]))
                
            self.G, self.has_groupvec = self._GenerateGroupMatrix(self.cids)
        
            # save everything to the database
            self.db.SaveSparseNumpyMatrix(self.STOICHIOMETRIC_TABLE_NAME, self.S)
            self.db.SaveSparseNumpyMatrix(self.GROUP_TABLE_NAME, self.G)
            self.db.SaveNumpyMatrix(self.GIBBS_ENERGY_TABLE_NAME, self.b.T)
            self.db.SaveNumpyMatrix(self.ANCHORED_TABLE_NAME, self.anchored.T)
            self.db.CreateTable(self.COMPOUND_TABLE_NAME, 'cid INT, name TEXT')
            for cid in self.cids:
                self.db.Insert(self.COMPOUND_TABLE_NAME, [cid, self.kegg.cid2name(cid)])
            self.db.CreateTable(self.UNIQUE_OBSERVATION_TABLE_NAME,
                                'row INT, id TEXT, type TEXT, url TEXT')
            for i in xrange(len(self.obs_ids)):
                self.db.Insert(self.UNIQUE_OBSERVATION_TABLE_NAME,
                               [i, self.obs_ids[i], self.obs_types[i], self.obs_urls[i]])
            self.db.Commit()

    def GetChemicalReactionEnergies(self, S, cids):
        return self._GetChemicalReactionEnergies(self.S.copy(), self.cids,
                                                 self.b.copy(), self.anchored,
                                                 S, cids)

    def _GetChemicalReactionEnergies(self, obs_S, obs_cids, obs_b, obs_anchored,
                                     est_S, est_cids):

        assert obs_S.shape[0] == len(obs_cids)
        assert obs_S.shape[1] == obs_b.shape[1]
        assert obs_S.shape[1] == obs_anchored.shape[1]
        assert est_S.shape[0] == len(est_cids)
        
        # Augment the observed stoichiometric matrix and the required S matrix
        # so that they will have the same number of rows and in the same order
        # of CIDs. Note that we must take the union of cids from both original
        # matrices.
        new_cids = list(set(est_cids).difference(obs_cids))
        all_cids = self.cids + new_cids
        obs_S = np.vstack([obs_S, np.matrix(np.zeros((len(new_cids), obs_S.shape[1])))])
        contributions, g_pgc, P_L_bad, P_L_pgc, G_resid = self._GetContributionData(
                                        obs_S, all_cids, obs_b, obs_anchored)
        
        new_est_S = np.matrix(np.zeros((len(all_cids), est_S.shape[1])))
        for c in abs(est_S).sum(1).nonzero()[0].flat:
            new_c = all_cids.index(est_cids[c])
            new_est_S[new_c, :] = est_S[c, :]
        est_S = new_est_S
        
        dG0_r = np.matrix(np.zeros((len(contributions), est_S.shape[1])))
        # now calculate the estimated dG0:
        for i, rowdict in enumerate(contributions):
            S_part = rowdict['P'] * est_S
            dG0_part = rowdict['g'] * S_part
            dG0_r[i, :] = dG0_part
            for j in xrange(est_S.shape[1]):
                r_str = UnifiedGroupContribution.row2string(S_part[:, j].round(10), all_cids)
                logging.debug("%s : %s, dG0 = %.1f" % (rowdict['name'], r_str, dG0_part))

        resid_S = P_L_bad * est_S
        resid_GS = P_L_pgc * G_resid * est_S

        resid = abs(resid_S).sum(0) + abs(resid_GS).sum(0)
        i_resid = list(np.where(resid > self.epsilon)[1].flat)
        dG0_r[:, i_resid] = np.nan
        return dG0_r

    def _GetContributionData(self, obs_S, obs_cids, obs_b, obs_anchored):
        assert obs_S.shape[0] == len(obs_cids)
        assert obs_S.shape[1] == obs_b.shape[1]
        assert obs_S.shape[1] == obs_anchored.shape[1]

        # (1)
        # use the anchored reactions to directly estimate the part of est_S
        # which is in their column-span, and normalize that part out from all matrices.
        anchored_cols = list(obs_anchored.nonzero()[1].flat)
        g_anch, P_C_anch, P_L_anch = LinearRegression.LeastSquaresProjection(
                            obs_S[:, anchored_cols], obs_b[:, anchored_cols])
        obs_b -= g_anch * P_C_anch * obs_S # subtract the contribution of anchored reactions to obs_b
        obs_S = P_L_anch * obs_S # project obs_S on the residual space
        obs_S[np.where(abs(obs_S) <= self.epsilon)] = 0

        
        # (2)
        # calculate the reactant contributions from obs_S and obs_b, and use that
        # to estimate the part which is in the column-space of NIST.
        g_prc, P_C_prc, P_L_prc = LinearRegression.LeastSquaresProjection(
                                                                obs_S, obs_b)

        # (3)
        # calculate the group contributions from obs_S and obs_b. Note that 
        # some reaction involve compounds that don't have groupvectors, and 
        # therefore are discarded from this step.
        G, has_groupvec = self._GenerateGroupMatrix(obs_cids)
        bad_compounds = list(np.where(has_groupvec == False)[0].flat)
        reactions_with_groupvec = []
        for i in xrange(obs_S.shape[1]):
            if np.all(abs(obs_S[bad_compounds, i]) < self.epsilon):
                reactions_with_groupvec.append(i)
        obs_GS = G.T * obs_S[:, reactions_with_groupvec]
        g_pgc, P_C_pgc, P_L_pgc = LinearRegression.LeastSquaresProjection(
                                    obs_GS, obs_b[:, reactions_with_groupvec])

        # calculate the total contributions
        contributions = []
        contributions.append({'g':g_anch, 'P':P_C_anch, 'name':'anchors'})
        contributions.append({'g':g_prc, 'P':P_C_prc * P_L_anch, 'name':'reactants'})
        contributions.append({'g':g_pgc * P_C_pgc * G.T, 'P':P_L_prc * P_L_anch, 'name':'groups'})
        
        #g_tot =  g_anch.copy()                              # anchored reaction contributions
        #g_tot += g_prc * P_C_prc * P_L_anch                 # reactant contributions
        #g_tot += g_pgc * P_C_pgc * G.T * P_L_prc * P_L_anch # group contributions
        
        P_L_bad = (P_L_prc * P_L_anch)[bad_compounds, :]    # projection for checking of we rely on compounds that have no groupvector
        G_resid = G.T * P_L_prc * P_L_anch                  # projection on the groupvector space

        #return g_tot, g_pgc, P_L_bad, P_L_pgc, G_resid
        return contributions, g_pgc, P_L_bad, P_L_pgc, G_resid
    
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
                      if S_row[c])
        r = Reaction("", sparse)
        return r.FullReactionString(show_cids=False)

#   def Estimate(self, S, b, r, cids, verbose=False):
#        """
#            Given
#                - group matrix: G (m x g)
#                - stoichiometric matrix: S (m x n)
#                - Gibbs energies: b (1 x n) 
#                - a new reaction: r (m x k)
#            Return:
#                The estimated Gibbs energy of 'r'
#        """
#        est = np.matrix(np.zeros((3, r.shape[1]))) * np.nan
#        
#        bad_compounds = np.where(self.has_groupvec == False)[0]
#        reactions_with_groupvec = []
#        for i in xrange(S.shape[1]):
#            if np.all(abs(S[bad_compounds, i]) < self.epsilon):
#                reactions_with_groupvec.append(i)
#        GS = self.G_used.T * S[:, reactions_with_groupvec]
#
#        try:
#            # calculate the contributions of compounds
#            g_S, PC_S, PL_S = LinearRegression.LeastSquaresProjection(S, b)
#            g_GS, PC_GS, PL_GS = LinearRegression.LeastSquaresProjection(GS, b[:, reactions_with_groupvec])
#        except np.linalg.linalg.LinAlgError:
#            return est
#
#        r_C = PC_S * r # the part of 'r' which is in the column-space of S
#        r_L = PL_S * r # the part of 'r' which is in the left-null-space of S
#
#
#        r_groupvec = self.G_used.T * r     # a group-vector representation of r
#        r_groupvec_L = self.G_used.T * r_L # a group-vector representation of r_L
#
#        est[0, :] = g_S * r_C
#        est[1, :] = g_GS * r_groupvec
#        est[2, :] = g_GS * r_groupvec_L + est[0, :]
#        
#        for k in xrange(r.shape[1]):
#            if verbose:
#                print "r =", UnifiedGroupContribution.row2string(r[:, k], self.cids)
#                print "r_C =", UnifiedGroupContribution.row2string(r_C[:, k], self.cids)
#                print "r_L =", UnifiedGroupContribution.row2string(r_L[:, k], self.cids)
#
#            if np.any(abs(r_L[:, k]) > self.epsilon):
#                # this reaction is not in the column space of S 
#                est[0, k] = np.nan
#
#            if (abs(r[bad_compounds, k]) > self.epsilon).any():
#                # this reaction involves compounds which don't have groupvectors
#                est[1, k] = np.nan
#            elif (abs(PL_GS * r_groupvec[:, k]) > self.epsilon).any():
#                # this reaction's groupvector is not in the column space of GS
#                est[1, k] = np.nan
#
#            if (abs(r_L[bad_compounds, k]) > self.epsilon).any():
#                # this reaction involves compounds which don't have groupvectors
#                est[2, k] = np.nan
#            elif (abs(PL_GS * r_groupvec_L[:, k]) > self.epsilon).any():
#                # this reaction's groupvector is not in the column space of GS
#                est[2, k] = np.nan
#        
#        return est
#    
#    def SqueezeData(self, normalize_anchors=True):
#        if normalize_anchors:
#            # now remove anchored data from S and leave only the data which will be 
#            # used for calculating the group contributions
#            anchored_cols = list(np.where(self.anchored==1)[1].flat)
#    
#            g, P_C, P_L = LinearRegression.LeastSquaresProjection(self.S[:, anchored_cols],
#                                                                  self.b[:, anchored_cols])
#    
#            # calculate the matrix and observations which are explained
#            # by the anchored reactions
#            self.S_anchored = P_C * self.S
#            self.b_anchored = g * self.S_anchored
#            self.cids_anchored = list(self.cids)
#            
#            # calculate the matrix and observations which are in the null-space
#            # of the anchored reactions. in other words, b_L are the residuals
#            # of the linear regression 
#            self.S -= self.S_anchored
#            self.b -= self.b_anchored
#            
#            # set epsilon-small values to absolute 0
#            self.S[np.where(abs(self.S) < self.epsilon)] = 0
#        else:
#            self.S_anchored = np.matrix(np.zeros(self.S.shape))
#            self.b_anchored = np.matrix(np.zeros(self.b.shape))
#            self.cids_anchored = list(self.cids)
#        
#        # removed zero rows (compounds) from S
#        used_cid_indices = list(np.nonzero(np.sum(abs(self.S), 1))[0].flat)
#        self.S_used = self.S[used_cid_indices, :]
#        self.G_used = self.G[used_cid_indices, :]
#        self.cids_used = [self.cids[i] for i in used_cid_indices]
#        self.has_groupvec_used = self.has_groupvec[used_cid_indices, :]
#
#    def Report(self, est):
#        if est.shape[0] == 1:
#            legend = ['UGC']
#        elif est.shape[0] == 3:
#            legend = ['PRC', 'PGC', 'UGC']
#        else:
#            raise ValueError()
#
#        resid = est - np.repeat(self.b.T, est.shape[0], 1).T
#        
#        fig = plt.figure(figsize=(6,6))
#        rms_list = []
#        for j, label in enumerate(legend):
#            used_reactions = list(np.where(np.isfinite(resid[j, :]))[1].flat)
#            x = list(self.b[0, used_reactions].flat)
#            y = list(resid[j, used_reactions].flat)
#            rms = rms_flat(y)
#            rms_list.append("RMSE(%s) = %.1f, N = %d" % (label, rms, len(x)))
#            plt.plot(x, y, '.', figure=fig, label=label)
#        
#        self.html_writer.write_ul(rms_list)
#        
#        plt.xlabel(r"$\Delta_r G^{'\circ}$ observed [kJ/mol]")
#        plt.ylabel(r"$\Delta_r G^{'\circ}$ residual [kJ/mol]")
#        plt.legend()
#        self.html_writer.embed_matplotlib_figure(fig)
#        
#        self.html_writer.insert_toggle(start_here=True, label="Show table")
#        rowdicts = []
#        for i in xrange(est.shape[1]):
#            rowdict = {}
#            rowdict['ID'] = self.obs_ids[i]
#            rowdict['type'] = self.obs_types[i]
#            rowdict['link'] = self.obs_urls[i]
#            #rowdict['reaction (anch)'] = UnifiedGroupContribution.row2hypertext(self.S_anchored[:, i], self.cids_anchored)
#            #rowdict['observed (anch)'] = "%.1f" % self.b_anchored[0, i]
#            rowdict['reaction'] = UnifiedGroupContribution.row2hypertext(self.S[:, i], self.cids)
#            rowdict['observed'] = "%.1f" % self.b[0, i]
#            for j, label in enumerate(legend):
#                rowdict[label] = "%.1f" % est[j, i]
#            rowdict['key'] = abs(resid[0, i])
#            rowdicts.append(rowdict)
#        #rowdicts.sort(key=lambda x:x['key'], reverse=True)
#        self.html_writer.write_table(rowdicts,
#            headers=['#', 'ID', 'type', 'reaction', 'observed'] + legend)
#            #headers=['#', 'reaction (anch)', 'observed (anch)', 'reaction', 'observed'] + legend)
#        self.html_writer.div_end()
#        
#    def Loo(self):
#        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')
#
#        n = self.S.shape[1]
#        resid = np.matrix(np.zeros((3, n))) * np.nan
#        for i in xrange(n):
#            if self.obs_types[i] == 'formation':
#                continue
#            if np.all(abs(self.S[:, i]) < self.epsilon): # empty reaction
#                continue
#            no_i = range(0, i) + range(i+1, n)
#            est = self.Estimate(self.S[:, no_i], self.b[:, no_i], 
#                                self.S[:, i:i+1], verbose=True)
#            resid[:, i:i+1] = est - self.b[0, i]
#            print 'dG0\' = %7.1f' % self.b[0, i],
#            print '| PRC = %5.1f' % resid[0, i],
#            print '| PGC = %5.1f' % resid[1, i],
#            print '| UGC = %5.1f' % resid[2, i],
#            print '| %s' % UnifiedGroupContribution.row2string(self.S[:, i], self.cids)        
#        self.Report(resid)
#        
#    def Fit(self):
#        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
#        est = self.Estimate(self.S, self.b, self.S)
#        resid = est - np.repeat(self.b.T, est.shape[0], 1).T
#        self.Report(resid)
    
    def Report(self, est):
        est = est.sum(0)
        finite = np.isfinite(est)
        resid = abs(self.b[finite] - est[finite])
        fig = plt.figure(figsize=(5,5), dpi=90)
        cdf(list(resid.flat), figure=fig)
        #plt.plot(self.b[finite].T, est[finite].T, '.', figure=fig)
        plt.title("RMSE = %.1f, N = %d" % (rms_flat(resid.flat), resid.shape[1]))
        plt.xlabel(r"$|\Delta_r G^{'\circ} obs - \Delta_r G^{'\circ} est|$ [kJ/mol]")
        plt.ylabel(r"CDF")
        self.html_writer.embed_matplotlib_figure(fig)

        rowdicts = []
        for i in xrange(self.b.shape[1]):
            rowdict = {}
            rowdict['row'] = i
            rowdict['reaction'] = UnifiedGroupContribution.row2hypertext(self.S[:, i], self.cids)
            rowdict['anchored'] = self.anchored[0, i]
            rowdict['obs'] = self.b[0, i]
            rowdict['est'] = est[0, i]
            if np.isfinite(est[0, i]):
                rowdict['|err|'] = abs(self.b[0, i] - est[0, i])
            else:
                rowdict['|err|'] = 0 
            rowdicts.append(rowdict)

        rowdicts.sort(key=lambda x:x['|err|'], reverse=True)            
        self.html_writer.insert_toggle(start_here=True, label="Show table")
        self.html_writer.write_table(rowdicts,
            headers=['row', 'reaction', 'anchored', 'obs', 'est', '|err|'], decimal=1)
        self.html_writer.div_end()
    
    def Fit(self):
        self.html_writer.write('<h2>Linear Regression Fit Analysis</h2>\n')
        est = self.GetChemicalReactionEnergies(self.S, self.cids)
        self.Report(est)

    def Loo(self):
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')
        n = self.S.shape[1]
        est = np.matrix(np.zeros((3, n))) * np.nan
        for i in xrange(n):
            no_i = range(0, i) + range(i+1, n)
            obs_S = self.S[:, no_i].copy()
            obs_anchored = self.anchored[0, no_i]
            obs_b = self.b[:, no_i].copy()
            est_S = self.S[:, i].copy()
            if np.all(abs(est_S) < self.epsilon): # empty reaction
                continue
            est[:, i] = self._GetChemicalReactionEnergies(obs_S, self.cids, obs_b, obs_anchored,
                                                          est_S, self.cids)
            logging.debug(UnifiedGroupContribution.row2string(self.S[:, i], self.cids))
            est_str = ' + '.join(['%.1f' % est[j, i] for j in xrange(est.shape[0])])
            logging.debug("obs = %.1f, est = %s" % (self.b[0, i], est_str))
        
        self.Report(est)

    def EstimateKeggCids(self):
        contributions, g_pgc, P_L_bad, P_L_pgc, G_resid = self._GetContributionData(
                        self.S.copy(), self.cids, self.b.copy(), self.anchored)
        
        g_tot = np.matrix(np.zeros((1, len(self.cids))))
        for rowdict in contributions:
            g_tot += rowdict['g'] * rowdict['P']
        
        diss = DissociationConstants.FromPublicDB()
        all_cids = sorted(self.kegg.get_all_cids())
        n_bad = P_L_bad.shape[0]
        n_pgc = P_L_pgc.shape[0]
        self.P_L_tot = np.matrix(np.zeros((n_bad + n_pgc, len(all_cids))))

        for c, cid in enumerate(all_cids):
            if cid not in self.cid2nH_nMg:
                self.cid2error[cid] = "No pKa data"
                continue
            
            nH, nMg = self.cid2nH_nMg[cid]
            if cid in self.cids:
                i = self.cids.index(cid)
                dG0 = g_tot[0, i]
                self.cid2source_string[cid] = "Unified Group Contribution"
                self.P_L_tot[0:n_bad, c] = P_L_bad[:, i] 
                self.P_L_tot[n_bad:, c] = P_L_pgc * G_resid[:, i]
            elif self.cid2groupvec[cid] is not None:
                gv = np.matrix(self.cid2groupvec[cid].Flatten())
                dG0 = float(g_pgc * gv.T)
                self.cid2source_string[cid] = "Group Contribution"
                self.P_L_tot[n_bad:, c] = P_L_pgc * gv.T
            else:
                self.cid2error[cid] = "no groupvector"
                continue
            
            diss_table = diss.GetDissociationTable(cid)
            if diss_table is not None:
                diss_table.SetFormationEnergyByNumHydrogens(
                    dG0=dG0, nH=nH, nMg=nMg)
                pmap = diss_table.GetPseudoisomerMap()
            else:
                pmap = PseudoisomerMap()
                pmap.Add(nH=nH, z=0, nMg=nMg, dG0=dG0)
            self.SetPseudoisomerMap(cid, pmap)
        
        conservation_rows = []        
        for i in xrange(self.P_L_tot.shape[0]):
            row_i = self.P_L_tot[i, :]
            c_active = sorted((abs(row_i) > self.epsilon).nonzero()[1].flat)
            if len(c_active) == 0:
                continue
            row_i = row_i * (1.0 / row_i[0, c_active[0]])
            row_i = row_i.round(10)
            
            # normalize reaction such that the coefficient of the smallest CID is 1
            sparse = dict((all_cids[c], row_i[0, c]) for c in c_active)
            if len(sparse) > 0:
                json_str = json.dumps(sparse)
                if i < n_bad:
                    conservation_rows.append(('missing structures and unknown reactant combination', json_str))
                else:
                    conservation_rows.append(('unknown reactant and group combination', json_str))
        
        self.db.CreateTable(self.CONSERVATIONS_TABLE_NAME, 'msg TEXT, json TEXT')
        conservation_rows = sorted(set(conservation_rows))
        for msg, json_str in conservation_rows:
            self.db.Insert(self.CONSERVATIONS_TABLE_NAME, [msg, json_str])
        self.ToDatabase(self.db, self.THERMODYNAMICS_TABLE_NAME,
                        self.ERRORS_TABLE_NAME)
        self.db.Commit()
    
    def init(self):
        if self.db.DoesTableExist(self.THERMODYNAMICS_TABLE_NAME):
            reader = self.db.DictReader(self.THERMODYNAMICS_TABLE_NAME)
            PsuedoisomerTableThermodynamics._FromDictReader(
                reader, self, label=None, name="Unified Group Contribution",
                warn_for_conflicting_refs=False)
    
            conservation_rows = []        
            for row in self.db.DictReader(self.CONSERVATIONS_TABLE_NAME):
                sparse = json.loads(row['json'])
                msg = row['msg']
                conservation_rows.append((msg, sparse))
    
            all_cids = sorted(self.kegg.get_all_cids())
            self.P_L_tot = np.matrix(np.zeros((len(conservation_rows), len(all_cids))))
            for i, (msg, sparse) in enumerate(conservation_rows):
                for cid, coeff in sparse.iteritems():
                    self.P_L_tot[i, all_cids.index(int(cid))] = float(coeff)
        else:
            self.LoadGroups(True)
            self.LoadObservations(True)
            self.LoadGroupVectors(True)
            self.LoadData(True)
            self.EstimateKeggCids()
    
    def GetTransfromedReactionEnergies(self, S, cids, pH=None, I=None, pMg=None, T=None):
        dG0_r = PsuedoisomerTableThermodynamics.GetTransfromedReactionEnergies(
                                                  self, S, cids, pH, I, pMg, T)

        # test to see if any of the reactions in S violate any conservation laws
        all_cids = sorted(self.kegg.get_all_cids())
        S_expanded = np.matrix(np.zeros((len(all_cids), S.shape[1])))
        for c, cid in enumerate(cids):
            i = all_cids.index(cid)
            S_expanded[i, :] = S[c, :]

        violations = abs(self.P_L_tot * S_expanded).sum(0) > self.epsilon
        dG0_r[violations] = np.nan        
        return dG0_r

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-d", "--from_database", action="store_true",
                          dest="from_database", default=False,
                          help="A flag for loading the data from the DB instead of "
                               "the CSV files (saves time but no debug information)")
    opt_parser.add_option("-t", "--train", action="store_true",
                          dest="train_only", default=False,
                          help="A flag for running the TRAIN only (without TEST)")
    return opt_parser
    
if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)

    options, _ = MakeOpts().parse_args(sys.argv)
    util._mkdir('../res')
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    html_writer = HtmlWriter('../res/ugc.html')
    
    ugc = UnifiedGroupContribution(db, html_writer)
    if options.train_only:
        ugc.LoadGroups(options.from_database)
        ugc.LoadObservations(options.from_database)
        ugc.LoadGroupVectors(options.from_database)
        ugc.LoadData(options.from_database)
        ugc.EstimateKeggCids()
        sys.exit(0)

    ugc.LoadGroups(True)
    ugc.LoadObservations(True)
    ugc.LoadGroupVectors(True)
    ugc.LoadData(True)
    ugc.init()

    r_list = []
    r_list += [Reaction.FromFormula("C00002 + C00001 = C00008 + C00009")]
    r_list += [Reaction.FromFormula("C00036 + C00024 = C00022 + C00083")]
    r_list += [Reaction.FromFormula("C00036 + C00100 = C00022 + C04348")]
    r_list += [Reaction.FromFormula("C01013 + C00010 + C00002 = C05668 + C00020 + C00013")]
    r_list += [Reaction.FromFormula("C00091 + C00005 = C00232 + C00010 + C00006")]
    r_list += [Reaction.FromFormula("C00002 + C00493 = C00008 + C03175")]
    
    kegg = Kegg.getInstance()
    S, cids = kegg.reaction_list_to_S(r_list)

    dG0_prime = ugc.GetTransfromedReactionEnergies(S, cids, pH=7, I=0.1)
    
    ln_conc = np.matrix(np.ones((1, S.shape[0]))) * np.log(0.001)
    if 1 in cids:
        ln_conc[0, cids.index(1)] = 0 # H2O should have a concentration of 1
    RT = R * default_T
    dGc_prime = dG0_prime + RT * ln_conc * S
    for i in xrange(len(r_list)):
        dG0 = ugc.GetChemicalReactionEnergies(S[:, i], cids)
        print "dG0 = %8.1f, dG0' = %8.1f,  dGc' = %8.1f,  %s" % \
            (dG0.sum(0), dG0_prime[0, i], dGc_prime[0, i],
             r_list[i].FullReactionString(show_cids=False))
        
        print ' + '.join('%.1f' % d for d in dG0.flat)
        
    #ugc.Fit()
    #ugc.Loo()