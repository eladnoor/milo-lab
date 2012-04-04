import numpy as np
from matplotlib.mlab import rms_flat
import matplotlib.pyplot as plt
from collections import defaultdict

from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.kegg_observation import KeggObervationCollection
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from pygibbs.kegg import Kegg
from toolbox.linear_regression import LinearRegression
import sys
from pygibbs.kegg_reaction import Reaction
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics,\
    AddConcentrationsToReactionEnergies
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
        
        self.FORMATION_ENERGY_FILENAME = '../data/thermodynamics/formation_energies.csv'

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
            self.obs_collection = KeggObervationCollection.FromDatabase(
                db=self.db,
                table_name=self.OBSERVATION_TABLE_NAME,
                transformed=self.transformed)
        else:
            logging.info("Reading observations from files")
            dissociation = self.GetDissociationConstants()
            self.obs_collection = KeggObervationCollection.FromFiles(
                html_writer=self.html_writer, 
                dissociation=dissociation,
                transformed=self.transformed,
                formation_energy_fname=self.FORMATION_ENERGY_FILENAME)
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
            # Here we simply use the dictionary self.cid2nH_nMg from KeggObervationCollection
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
            logging.info("Reading group matrices from database")
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
            for rowdict in self.db.DictReader(self.UNIQUE_OBSERVATION_TABLE_NAME):
                self.obs_ids.append(rowdict['id'])
                self.obs_types.append(rowdict['type'])
                self.obs_urls.append(rowdict['url'])
        else:
            logging.info("Calculating group matrices")
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
                self.obs_types.append(', '.join(set([obs.obs_type for obs in obs_list])))
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

    def _GetContributionData(self, obs_S, obs_cids, obs_b, obs_anchored):
        assert obs_S.shape[0] == len(obs_cids)
        assert obs_S.shape[1] == obs_b.shape[1]
        assert obs_S.shape[1] == obs_anchored.shape[1]

        # (1)
        # use the anchored reactions to directly estimate the part of est_S
        # which is in their column-span, and normalize that part out from all matrices.
        anchored_cols = list(obs_anchored.nonzero()[1].flat)
        if anchored_cols:
            g_anch, P_C_anch, P_L_anch = LinearRegression.LeastSquaresProjection(
                                obs_S[:, anchored_cols], obs_b[:, anchored_cols])
            obs_b -= g_anch * P_C_anch * obs_S # subtract the contribution of anchored reactions to obs_b
            obs_S = P_L_anch * obs_S # project obs_S on the residual space
        else:
            g_anch = np.matrix(np.zeros((1, obs_S.shape[0])))
            P_C_anch = np.matrix(np.zeros((obs_S.shape[0], obs_S.shape[0])))
            P_L_anch = np.matrix(np.eye(obs_S.shape[0]))
        
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
        result_dict = {}
        result_dict['names'] = ['anchors', 'reactants', 'groups']
        result_dict['contributions'] = [g_anch, g_prc, g_pgc * P_C_pgc * G.T]
        result_dict['group_contributions'] = g_pgc
        result_dict['column_spaces'] = [P_C_anch, P_C_prc, P_C_pgc]
        result_dict['null_spaces'] = [P_L_anch, P_L_prc, P_L_pgc]
        result_dict['projections'] = [P_C_anch,
                                      P_C_prc * P_L_anch,
                                      P_L_prc * P_L_anch]
        
        result_dict['total_contributions'] = np.matrix(np.zeros((1, len(obs_cids))))
        for g, S in zip(result_dict['contributions'], result_dict['projections']):
            result_dict['total_contributions'] += g * S
        
        # conservation laws that check if we rely on compounds that have no groupvector
        P_L_bad = (P_L_prc * P_L_anch)[bad_compounds, :]
        
        # projection of reactions to the residual groupvector space
        G_resid = P_L_prc * P_L_anch * G

        result_dict['bad_conservations'] = P_L_bad 
        result_dict['pgc_conservations'] = P_L_pgc
        result_dict['pgc_groupvectors'] = G_resid
        result_dict['conservations'] = np.vstack([P_L_bad, (G_resid * P_L_pgc).T])

        return result_dict

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
        result_dict = self._GetContributionData(obs_S, all_cids,
                                                     obs_b, obs_anchored)
        
        new_est_S = np.matrix(np.zeros((len(all_cids), est_S.shape[1])))
        for c in abs(est_S).sum(1).nonzero()[0].flat:
            new_c = all_cids.index(est_cids[c])
            new_est_S[new_c, :] = est_S[c, :]
        est_S = new_est_S
        
        n_cont = len(result_dict['names'])
        dG0_r = np.matrix(np.zeros((n_cont, est_S.shape[1])))
        
        # each value i, j in the matrix 'parts' is the norm2 of the
        # projection of the reaction S[:, j] on the subspace corresponding
        # to result_dict['projections'][i]
        parts = np.matrix(np.zeros((n_cont+1, est_S.shape[1])))
        
        # now calculate the estimated dG0:
        for i in xrange(n_cont):
            S_part = result_dict['projections'][i] * est_S
            parts[i, :] = np.sqrt(np.diag(S_part.T * S_part)) # the norm2 of each column in S_part
            dG0_part = result_dict['contributions'][i] * S_part
            dG0_r[i, :] = dG0_part

            # only print the debug info if there is one target reaction
            # otherwise, there is too much info to printout
            for j in xrange(est_S.shape[1]):
                r_str = UnifiedGroupContribution.row2string(S_part[:, j].round(10), all_cids)
                logging.debug("%s : %s, dG0 = %.1f (%.2g)" %
                              (result_dict['names'][i], r_str, dG0_r[i, j], parts[i, j]))

        parts[-1, :] = abs(result_dict['conservations'] * est_S).sum(0)
        i_resid = list((parts[-1, :] > self.epsilon).nonzero()[1].flat)
        dG0_r[:, i_resid] = np.nan
        
        # calculate also dG0_r for the PGC method without using PRC (i.e.
        # classic group contribution).
        dG0_r_pgc = dG0_r[0, :].copy() # the anchored part is the same as before
        S_part_pgc = result_dict['null_spaces'][0] * est_S
        dG0_r_pgc += result_dict['contributions'][2] * S_part_pgc
        resid_pgc = abs(result_dict['null_spaces'][2] * result_dict['pgc_groupvectors'].T * S_part_pgc).sum(0)
        i_resid_pgc = list((resid_pgc > self.epsilon).nonzero()[1].flat)
        dG0_r_pgc[0, i_resid_pgc] = np.nan
        
        return dG0_r, parts, dG0_r_pgc
    
    def EstimateKeggCids(self):
        result_dict = self._GetContributionData(self.S.copy(), self.cids,
                                                self.b.copy(), self.anchored)
        
        g_pgc = result_dict['group_contributions']
        P_L_bad = result_dict['bad_conservations']
        P_L_pgc = result_dict['pgc_conservations']
        G_resid = result_dict['pgc_groupvectors']
        g_tot = result_dict['total_contributions']
        
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
                self.P_L_tot[:n_bad, c] = P_L_bad[:, i] 
                self.P_L_tot[n_bad:, c] = P_L_pgc * G_resid[i, :].T
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

    def Report(self, est, title):
        self.html_writer.write('</br><b>%s</b><br>\n' % title)

        finite = np.isfinite(est)
        resid = abs(self.b[finite] - est[finite])
        fig = plt.figure(figsize=(5,5), dpi=60)
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
            rowdict['type'] = self.obs_types[i]
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
            headers=['row', 'type', 'reaction', 'anchored', 'obs', 'est', '|err|'], decimal=1)
        self.html_writer.div_end()
    
    def Fit(self):
        dG0_r_ugc, _, dG0_r_pgc = self.GetChemicalReactionEnergies(self.S, self.cids)
        self.Report(dG0_r_ugc.sum(0), 'UGC - regression fit')
        self.Report(dG0_r_pgc, 'PGC - regression fit')

    def Loo(self):
        n = self.S.shape[1]
        dG0_r_ugc = np.matrix(np.zeros((3, n))) * np.nan
        dG0_r_pgc = np.matrix(np.zeros((1, n))) * np.nan

        rowdicts = []
        class2ugc_err = defaultdict(list)
        class2pgc_err = defaultdict(list)
        for i in xrange(n):
            if self.obs_types[i] != 'reaction':
                continue
            if self.anchored[0, i]:
                continue
            if abs(self.S[:, i]).sum(0) < self.epsilon: # empty reaction
                continue

            no_i = range(0, i) + range(i+1, n)
            obs_S = self.S[:, no_i].copy()
            obs_anchored = self.anchored[0, no_i]
            obs_b = self.b[:, no_i].copy()
            est_S = self.S[:, i].copy()
            dG0_r_ugc[:, i], parts, dG0_r_pgc[0, i] = self._GetChemicalReactionEnergies(
                obs_S, self.cids, obs_b, obs_anchored, est_S, self.cids)

            if parts[3, 0] > self.epsilon:
                classification = 'kernel'
            elif parts[1, 0] > self.epsilon and parts[2, 0] > self.epsilon:
                classification = 'PRC + PGC'
            elif parts[1, 0] > self.epsilon:
                classification = 'PRC'
            elif parts[2, 0] > self.epsilon:
                classification = 'PGC'
            else:
                classification = 'anchored'
            
            est_b = float(dG0_r_ugc[:, i].sum(0))
            ugc_err = self.b[0, i] - est_b
            class2ugc_err[classification].append(ugc_err)
            
            pgc_err = self.b[0, i] - dG0_r_pgc[0, i]
            class2pgc_err[classification].append(pgc_err)

            rowdict = {}
            rowdict['row'] = i
            rowdict['type'] = self.obs_types[i]
            rowdict['reaction'] = UnifiedGroupContribution.row2hypertext(self.S[:, i], self.cids)
            rowdict['obs'] = self.b[0, i]
            rowdict['est'] = est_b
            rowdict['est(PGC)'] = dG0_r_pgc[0, i]
            if np.isfinite(ugc_err):
                rowdict['|err|'] = abs(ugc_err)
            else:
                rowdict['|err|'] = 0
            rowdict['est_ANCH'] = dG0_r_ugc[0, i]
            rowdict['est_PRC'] = dG0_r_ugc[1, i]
            rowdict['est_PGC'] = dG0_r_ugc[2, i]
            rowdict['part_ANCH'] = parts[0, 0]
            rowdict['part_PRC'] = parts[1, 0]
            rowdict['part_PGC'] = parts[2, 0]
            rowdict['part_NULL'] = parts[3, 0]
            rowdict['class'] = classification
            rowdicts.append(rowdict)
            
        class_errors = []
        for classification in class2ugc_err.keys():
            ugc_err_list = class2ugc_err[classification]
            pgc_err_list = class2pgc_err[classification]
            class_errors.append('%s: N = %d, rmse(UGC) = %.1f kJ/mol, rmse(PGC) = %.1f kJ/mol' % 
                                (classification, len(ugc_err_list),
                                 rms_flat(ugc_err_list), rms_flat(pgc_err_list)))
        
        self.Report(dG0_r_ugc.sum(0), 'UGC - Leave one out')
        self.Report(dG0_r_pgc, 'PGC - Leave one out')

        rowdicts.sort(key=lambda x:x['|err|'], reverse=True)            
        self.html_writer.write('<h2>Linear Regression Leave-One-Out Analysis</h2>\n')
        self.html_writer.insert_toggle(start_here=True, label="Show table")
        self.html_writer.write_ul(class_errors)
        self.html_writer.write_table(rowdicts,
            headers=['row', 'type', 'reaction', 'class', 'obs', 'est',
                     'est(PGC)', '|err|', 'est_ANCH', 'est_PRC', 'est_PGC',
                     'part_ANCH', 'part_PRC', 'part_PGC', 'part_NULL'], decimal=1)
        self.html_writer.div_end()

    def init(self):
        if self.db.DoesTableExist(self.THERMODYNAMICS_TABLE_NAME):
            logging.info('Reading thermodynamic data from database')
            reader = self.db.DictReader(self.THERMODYNAMICS_TABLE_NAME)
            PsuedoisomerTableThermodynamics._FromDictReader(
                reader, self, label=None, name="Unified Group Contribution",
                warn_for_conflicting_refs=False)
    
            conservation_rows = []        
            for row in self.db.DictReader(self.CONSERVATIONS_TABLE_NAME):
                sparse = dict((int(cid), coeff) for (cid, coeff) in json.loads(row['json']).iteritems())
                msg = row['msg']
                conservation_rows.append((msg, sparse))
    
            logging.info('Reading conservation matrix data from database')
            all_cids = sorted(self.kegg.get_all_cids())
            cid_dict = dict((cid, i) for (i, cid) in enumerate(all_cids))
            self.P_L_tot = np.matrix(np.zeros((len(conservation_rows), len(all_cids))))
            for i, (msg, sparse) in enumerate(conservation_rows):
                for cid, coeff in sparse.iteritems():
                    self.P_L_tot[i, cid_dict[cid]] = float(coeff)
        else:
            self.LoadGroups(True)
            self.LoadObservations(True)
            self.LoadGroupVectors(True)
            self.LoadData(True)
            self.EstimateKeggCids()
    
    def GetTransfromedReactionEnergies(self, S, cids, pH=None, I=None, pMg=None, T=None, conc=1):
        dG0_r = PsuedoisomerTableThermodynamics.GetTransfromedReactionEnergies(
                                                  self, S, cids, pH, I, pMg, T)
        if conc != 1:
            pH, I, pMg, T = self.GetConditions(pH, I, pMg, T)
            dG0_r += AddConcentrationsToReactionEnergies(S, cids, T, conc)

        # test to see if any of the reactions in S violate any conservation laws
        all_cids = sorted(self.kegg.get_all_cids())
        S_expanded = np.matrix(np.zeros((len(all_cids), S.shape[1])))
        for c, cid in enumerate(cids):
            i = all_cids.index(cid)
            S_expanded[i, :] = S[c, :]

        violations = abs(self.P_L_tot * S_expanded).sum(0) > self.epsilon
        dG0_r[violations] = np.nan        
        return dG0_r
    
    def SaveDataToMatfile(self):
        np.savetxt(fname='../res/ugc_S.txt', X=self.S, fmt="%g", delimiter=',', newline='\n')
        np.savetxt(fname='../res/ugc_b.txt', X=self.b, fmt="%g", delimiter=',', newline='\n')
        np.savetxt(fname='../res/ugc_cids.txt', X=np.array(self.cids), fmt="%d", delimiter=',', newline='\n')
        np.savetxt(fname='../res/ugc_anchored.txt', X=self.anchored, fmt="%d", delimiter=',', newline='\n')
    
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-g", "--recalc_groups", action="store_true",
                          dest="recalc_groups", default=False,
                          help="A flag for loading the group definitions from the CSV"
                               " instead of the DB")
    opt_parser.add_option("-o", "--recalc_observations", action="store_true",
                          dest="recalc_observations", default=False,
                          help="A flag for loading the dG observations from the CSV"
                               " instead of the DB")
    opt_parser.add_option("-v", "--recalc_groupvectors", action="store_true",
                          dest="recalc_groupvectors", default=False,
                          help="A flag for recalculating the group vectors"
                               " instead of loading them from the DB")
    opt_parser.add_option("-m", "--recalc_matrices", action="store_true",
                          dest="recalc_matrices", default=False,
                          help="A flag for recalculating the group matrices"
                               " instead of loading them from the DB")
    opt_parser.add_option("-t", "--train", action="store_true",
                          dest="train", default=False,
                          help="A flag for running the TRAIN")
    opt_parser.add_option("-e", "--test", action="store_true",
                          dest="test", default=False,
                          help="A flag for running the TEST")
    opt_parser.add_option("-l", "--leave_one_out", action="store_true",
                          dest="loo", default=False,
                          help="A flag for running the Leave One Out analysis")
    return opt_parser
    
if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)

    options, _ = MakeOpts().parse_args(sys.argv)
    util._mkdir('../res')
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    html_writer = HtmlWriter('../res/ugc.html')
    
    ugc = UnifiedGroupContribution(db, html_writer)
    ugc.LoadGroups(FromDatabase=(not options.recalc_groups))
    ugc.LoadObservations(FromDatabase=(not options.recalc_observations))
    ugc.LoadGroupVectors(FromDatabase=(not options.recalc_groupvectors))
    ugc.LoadData(FromDatabase=(not options.recalc_matrices))
    ugc.SaveDataToMatfile()
    sys.exit(0)
    
    if options.train:
        ugc.EstimateKeggCids()
        sys.exit(0)
    else:
        ugc.init()

    if options.test:
        r_list = []
#        r_list += [Reaction.FromFormula("C00002 + C00001 = C00008 + C00009")]
#        r_list += [Reaction.FromFormula("C00036 + C00024 = C00022 + C00083")]
#        r_list += [Reaction.FromFormula("C00036 + C00100 = C00022 + C00683")]
#        r_list += [Reaction.FromFormula("C01013 + C00010 + C00002 = C05668 + C00020 + C00013")]
#        r_list += [Reaction.FromFormula("C00091 + C00005 = C00232 + C00010 + C00006")]
#        r_list += [Reaction.FromFormula("C00002 + C00493 = C00008 + C03175")]
#        r_list += [Reaction.FromFormula("C00243 + C00125 = C05403 + C00126")]
#        r_list += [Reaction.FromFormula("2 C00206 = C00360 + C00131")]
        r_list += [Reaction.FromFormula("C04171 + C00003 = C00196 + C00080 + C00004")]
        
        
        kegg = Kegg.getInstance()
        S, cids = kegg.reaction_list_to_S(r_list)
    
        dG0_prime = ugc.GetTransfromedReactionEnergies(S, cids, pH=7.0, I=0.15)
        
        ln_conc = np.matrix(np.ones((1, S.shape[0]))) * np.log(0.001)
        if 1 in cids:
            ln_conc[0, cids.index(1)] = 0 # H2O should have a concentration of 1
        RT = R * default_T
        dGc_prime = dG0_prime + RT * ln_conc * S
        for i in xrange(len(r_list)):
            r_list[i].Balance()
            dG0_r, parts, dG0_r_pgc = ugc.GetChemicalReactionEnergies(S[:, i], cids)
            print r_list[i].FullReactionString(show_cids=False)
            print ('UGC: dG0 = %.1f = ' % dG0_r.sum(0)) + ' + '.join('%.1f' % d for d in dG0_r.flat)
            print ('PGC: dG0 = %.1f = ' % dG0_r_pgc.sum(0)) + ' + '.join('%.1f' % d for d in dG0_r_pgc.flat)
            print "dG0' = %.1f, dGc' = %.1f" % \
                (dG0_prime[0, i], dGc_prime[0, i])
    
    if options.loo:
        ugc.Fit()
        ugc.Loo()