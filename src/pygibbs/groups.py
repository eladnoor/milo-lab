#!/usr/bin/python

import logging, types, json, sys
from collections import defaultdict

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rms_flat
from pygibbs.thermodynamic_constants import R, default_pH, default_T,\
    dG0_f_Mg, default_I, default_pMg, RedoxCarriers,\
    symbol_d_G0, symbol_d_G0_prime
from pygibbs.thermodynamics import MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics, AddConcentrationsToReactionEnergies
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.groups_data import Group, GroupsData
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from pygibbs.kegg_observation import KeggObervationCollection
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamic_errors import MissingReactionEnergy
from pygibbs.group_vector import GroupVector
from optparse import OptionParser
from toolbox import util

class GroupMissingTrainDataError(Exception):
    
    def __init__(self, value, msg, kernel_rows=[]):
        """
            value - the estimated value for this group vector
            msg   - the error message
            kernel_rows - the kernel rows which are not orthogonal to the input group vector
        """
        self.value = value
        self.msg = msg
        self.kernel_rows = kernel_rows
    
    def __str__(self):
        if type(self.msg) == types.StringType:
            return self.msg
        else:
            return repr(self.msg)
    
    def Explain(self, gc):
        missing_single_groups = []
        ker = abs(gc.group_nullspace) > 1e-10
        for i in self.kernel_rows:
            nonzero_columns = list(ker[:, i].nonzero()[0].flat)
            if len(nonzero_columns) == 1:
                missing_single_groups.append(nonzero_columns[0])
            else:
                return 'contains group combinations that are not covered by ' + \
                    'training data'
                
        return 'contains missing groups: ' + ", ".join(
            [str(gc.groups_data.all_groups[j]) for j in missing_single_groups])
        
class GroupContribution(PsuedoisomerTableThermodynamics):    

    def __init__(self, db, html_writer=None, transformed=False):
        """Construct a GroupContribution instance.
        
        Args:
            db: the database handle to read from.
            html_writer: the HtmlWriter to write to.
            kegg: a Kegg instance if you don't want to use the default one.
        """
        PsuedoisomerTableThermodynamics.__init__(self, name="Group Contribution")
        self.db = db
        self.html_writer = html_writer or NullHtmlWriter()
        self.dissociation = None
        self.transformed = transformed
        
        self.epsilon = 1e-10

        self.kegg = Kegg.getInstance()
        self.bounds = deepcopy(self.kegg.cid2bounds)

        self.group_nullspace = None
        self.group_contributions = None
        self.obs_collection = None
        
        self.cid2error = {}
        self.cid2groupvec = None

        if transformed:
            prefix = 'bgc'
        else:
            prefix = 'pgc'
        
        self.OBSERVATION_TABLE_NAME = prefix + '_observations'
        self.GROUPVEC_TABLE_NAME = prefix + '_groupvector'
        self.NULLSPACE_TABLE_NAME = prefix + '_nullspace'
        self.CONTRIBUTION_TABLE_NAME = prefix + '_contribution'
        self.REGRESSION_TABLE_NAME = prefix + '_regression'
        
        self.THERMODYNAMICS_TABLE_NAME = prefix + '_pseudoisomers'
        self.STOICHIOMETRIC_MATRIX_TABLE_NAME = prefix + '_stoichiometry'
        self.ANCHORED_CONTRIBUTIONS_TALBE_NAME = prefix + '_anchored_g'
        self.ANCHORED_CIDS_TABLE_NAME = prefix + '_anchored_cids'
        self.ANCHORED_P_L_TALBE_NAME = prefix + '_anchored_P_L'
        
    def GetDissociationConstants(self):
        """
            Since loading the pKas takes time, this function is a lazy initialization
            of self.dissociation.
        """
        if self.dissociation is None:
            self.dissociation = DissociationConstants.FromPublicDB()
        return self.dissociation
    
    def GetDissociationTable(self, cid):
        return self.GetDissociationConstants().GetDissociationTable(cid,
                                                    create_if_missing=False)
    
    def WriteDataToJSON(self, json_fname, kegg):
        formations = []
        for row in self.db.DictReader('gc_cid2error'):
            h = {}
            h['cid'] = row['cid']
            try:
                h['name'] = kegg.cid2name(h['cid'])
            except KeyError:
                h['name'] = None
            try:
                h['inchi'] = kegg.cid2inchi(h['cid'])
            except KeyError:
                h['inchi'] = None
            h['num_electrons'] = kegg.cid2num_electrons(h['cid'])
            h['source'] = self.cid2source_string.get(row['cid'], None)
            h['species'] = []
            try:
                for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(row['cid']).ToMatrix():
                    h['species'].append({"nH":nH, "z":z, "nMg":nMg, "dG0_f":dG0})
            except MissingCompoundFormationEnergy:
                h['error'] = row['error']

            formations.append(h)

        json_file = open(json_fname, 'w')
        json_file.write(json.dumps(formations, indent=4))
        json_file.close()
        
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

            # When using non-transformed energies, it is very important for
            # the group vectors of each compound to represent the correct
            # pseudoisomer (same nH and nMg used in the reverse transform and
            # the list of formation energies). Here we use the dictionary 
            # self.cid2nH_nMg that is copied from GroupObervationCollection
            self.cid2nH_nMg = self.obs_collection.cid2nH_nMg
            if self.cid2nH_nMg is None:
                # Use I = 0 mM and pMg = 14 since these are the conditions used
                # to determine the most abundant pseudoisomer in DissociationConstants
                # (that is the only option in ChemAxon).
                # Since we rely on that table for the molecular structures,
                # we must be consistent with it here.
                self.cid2nH_nMg = dissociation.GetCid2nH_nMg(
                                            pH=self.pH, I=0, pMg=14, T=self.T)

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

    def Train(self):
        logging.info("Calculating the linear regression data")
        cids, S, b, anchored = self.obs_collection.GetStoichiometry()
        
        anchored_cols = list(np.where(anchored==1)[1].flat)
        # now remove anchored data from S and leave only the data which will be 
        # used for calculating the group contributions
        g, P_C, P_L = LinearRegression.LeastSquaresProjection(S[:, anchored_cols],
                                                              b[:, anchored_cols])
        self.anchored_cids = cids
        self.anchored_contributions = g * P_C
        self.anchored_P_L = P_L
        self.anchored_P_L[abs(self.anchored_P_L) <= self.epsilon] = 0

        b -= self.anchored_contributions * S
        S = self.anchored_P_L * S

        # set epsilon-small values to absolute 0
        S[np.where(abs(S) <= self.epsilon)] = 0
        
        # removed zero rows (compounds) from S
        used_cid_indices = set(np.nonzero(np.sum(abs(S), 1))[0].flat)
        for i_cid, cid in enumerate(cids):
            if self.cid2groupvec[cid] is None:
                used_cid_indices.difference_update([i_cid])
                for i_obs in np.nonzero(S[i_cid, :])[1].flat:
                    logging.warning("%s is removed because C%05d has no group vector, "
                                    "but is still part of the final stoichiometric matrix"
                                    % (self.obs_collection.observations[i_obs].obs_id,
                                       cid))
                    S[:, i_obs] = 0

        used_cid_indices = sorted(used_cid_indices)
        S = S[used_cid_indices, :]

        n_groups = len(self.groups_data.GetGroupNames()) # number of groups
        G = np.matrix(np.zeros((len(used_cid_indices), n_groups)))
        for i, i_cid in enumerate(used_cid_indices):
            G[i, :] = self.cid2groupvec[cids[i_cid]].Flatten()

        GS = G.T * S

        # 'unique' the rows GS. For each set of rows that is united,
        # the Y-value for the new row is the average of the corresponding Y-values.
        unique_GS, col_mapping = LinearRegression.ColumnUnique(GS, remove_zero=True)
        unique_b = np.matrix(np.zeros((1, unique_GS.shape[1])))
        unique_obs_types = []
        unique_obs_ids = []
        for i, old_indices in sorted(col_mapping.iteritems()):
            unique_b[0, i] = np.mean(b[0, old_indices])
            obs_list = [self.obs_collection.observations[j] for j in old_indices]
            unique_obs_types.append(obs_list[0].obs_type) # take the type of the first one (not perfect...)
            unique_obs_ids.append(', '.join([obs.obs_id for obs in obs_list]))            
        
        self.group_matrix = unique_GS
        self.obs_values = unique_b
        self.obs_ids = unique_obs_ids
        self.obs_types = unique_obs_types

        logging.info("Performing linear regression")
        self.group_contributions, self.group_nullspace = \
            LinearRegression.LeastSquares(self.group_matrix, self.obs_values)
        
        logging.info("Storing the group contribution data in the database")
        self.SaveContributionsToDB()

    def AnalyzeResiduals(self):
        GS = np.dot(self.G.T, self.S)
        # Write the analysis of residuals:
        # I am not sure if this analysis should be done before "uniquing"
        # the rows of S or after. The observation residual is much smaller
        # in the latter case, since intra-reaction noise is averaged.
        _P_R1, P_N1 = LinearRegression.RowProjection(self.S)
        _P_R2, P_N2 = LinearRegression.RowProjection(GS)
        
        r_obs = np.linalg.norm(np.dot(self.gibbs_values, P_N1))
        r_est = np.linalg.norm(np.dot(self.gibbs_values, P_N2 - P_N1))
        r_tot = np.linalg.norm(np.dot(self.gibbs_values, P_N2))
        
        self.html_writer.write('</br><b>Analysis of residuals:<b>\n')
        self.html_writer.insert_toggle(start_here=True)
        residual_text = ['r<sub>observation</sub> = %.2f kJ/mol' % r_obs,
                         'r<sub>estimation</sub> = %.2f kJ/mol' % r_est,
                         'r<sub>total</sub> = %.2f kJ/mol' % r_tot]
        self.html_writer.write_ul(residual_text)
        self.html_writer.div_end()

    def SaveContributionsToDB(self):
        self.db.CreateTable(self.ANCHORED_CIDS_TABLE_NAME, 'cid INT')
        for cid in self.anchored_cids:
            self.db.Insert(self.ANCHORED_CIDS_TABLE_NAME, [cid])
        self.db.SaveNumpyMatrix(self.ANCHORED_CONTRIBUTIONS_TALBE_NAME,
                                self.anchored_contributions.T)
        self.db.SaveSparseNumpyMatrix(self.ANCHORED_P_L_TALBE_NAME,
                                      self.anchored_P_L)
        
        # write a table of the group contributions
        self.db.CreateTable(self.CONTRIBUTION_TABLE_NAME,
                            'gid INT, name TEXT, protons INT, charge INT, '
                            'nMg INT, dG0_gr REAL')
        for j, group_name in enumerate(self.groups_data.GetGroupNames()):
            dG0_gr = self.group_contributions[0, j]
            
            if self.transformed:
                nH, z, nMg = None, None, None
            else:
                group = self.groups_data.all_groups[j]
                nH, z, nMg = group.hydrogens, group.charge, group.nMg

            self.db.Insert(self.CONTRIBUTION_TABLE_NAME,
                [j, group_name, nH, z, nMg, dG0_gr])
            
        self.db.SaveSparseNumpyMatrix(self.NULLSPACE_TABLE_NAME,
                                      self.group_nullspace)
        self.db.Commit()
            
    def LoadContributionsFromDB(self):
        logging.info("loading the group contribution data from the database")
        self.anchored_cids = []
        for row in self.db.DictReader(self.ANCHORED_CIDS_TABLE_NAME):
            self.anchored_cids.append(row['cid'])
        self.anchored_contributions = self.db.LoadNumpyMatrix(self.ANCHORED_CONTRIBUTIONS_TALBE_NAME).T
        self.anchored_P_L = self.db.LoadSparseNumpyMatrix(self.ANCHORED_P_L_TALBE_NAME)

        self.group_contributions = []
        for row in self.db.DictReader(self.CONTRIBUTION_TABLE_NAME):
            self.group_contributions.append(row['dG0_gr'])
        self.group_contributions = np.matrix([self.group_contributions])
        self.group_nullspace = self.db.LoadSparseNumpyMatrix(self.NULLSPACE_TABLE_NAME)
                        
    def does_table_exist(self, table_name):
        for unused_ in self.db.Execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def groupvec2val(self, groupvec):
        if self.group_contributions == None or self.group_nullspace == None:
            raise Exception("You need to first Train the system before using it to estimate values")

        gv = np.matrix(groupvec.Flatten()).T
        val = float(self.group_contributions * gv)
        v = self.group_nullspace.T * gv
        k_list = list(np.where(abs(v) > self.epsilon)[0].flat)
        if k_list:
            raise GroupMissingTrainDataError(val, "can't estimate because the input "
                "is not in the column-space of the group matrix", k_list)
        return val
    
    def WriteRegressionReport(self, T=default_T, pH=default_pH):        
        rowdicts = []
        for i in xrange(self.group_matrix.shape[1]):
            groupvec = GroupVector(self.groups_data, self.group_matrix[:, i])
            rowdict = {'#': i, 'ID': self.obs_ids[i]}
            rowdict[self.obs_collection.gibbs_symbol] = '%.1f' % self.obs_values[0, i]
            rowdict['Group Vector'] = str(groupvec)
            rowdicts.append(rowdict)

        self.html_writer.write('</br><b>Regression report</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.html_writer.write('<font size="1">\n')
        self.html_writer.write_ul(['observations: %d' % self.group_matrix.shape[1],
           'groups: %d' % self.group_matrix.shape[0],
           'rank: %d' % LinearRegression.MatrixRank(self.group_matrix)])
        self.html_writer.write_table(rowdicts, 
                        headers=['#', 'ID', 'Group Vector', self.obs_collection.gibbs_symbol])
        self.html_writer.write('</font>\n')
        self.html_writer.div_end()
        
        self.html_writer.write('</br><b>Group Contributions</b>\n')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)
        self.html_writer.write('</br><font size="1">\n')
        rowdicts = []
        if self.transformed:
            headers = ["#", "Group Name",
                self.obs_collection.gibbs_symbol, 
                "acid-base", "formation", "reaction"]
        else:
            headers = ["#", "Group Name",
                       "nH", "charge", "nMg",
                       self.obs_collection.gibbs_symbol, 
                       "acid-base", "formation", "reaction"]
        group_names = self.groups_data.GetGroupNames()
        for j, dG0_gr in enumerate(self.group_contributions.flat):
            obs_lists_dict = defaultdict(list)
            for k in self.group_matrix[j, :].nonzero()[1].flat:
                obs_lists_dict[self.obs_types[k]].append(self.obs_ids[k])
            d = {"#":"%d" % j, "Group Name":group_names[j], 
                 self.obs_collection.gibbs_symbol:"%.1f" % dG0_gr}
            for k, v in obs_lists_dict.iteritems():
                d[k] = ' | '.join(v)
            if not self.transformed:
                group = self.groups_data.all_groups[j]
                d["nH"] = group.hydrogens
                d["charge"] = group.charge
                d["nMg"] = group.nMg
            rowdicts.append(d)
        self.html_writer.write_table(rowdicts, headers)
        self.html_writer.write('</font>\n')
        self.html_writer.div_end()

    def AnalyzeTrainingSet(self):
        n_obs = self.group_matrix.shape[1]
        rowdicts = []
        fit_results = np.dot(self.group_contributions, self.group_matrix)
        residuals = fit_results - self.obs_values
        
        if self.transformed:
            sym = symbol_d_G0_prime
        else:
            sym = symbol_d_G0
        for i in xrange(n_obs):
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.obs_types[i] not in ['formation', 'reaction']:
                continue

            rowdict = {'Observation':self.obs_ids[i]}
            rowdict[sym + ' (obs)'] = self.obs_values[0, i]
            rowdict[sym + ' (fit)'] = fit_results[0, i]
            rowdict[sym + ' (res)'] = residuals[0, i]
            rowdict['LOO ' + sym + ' (fit)'] = np.nan
            rowdict['LOO ' + sym + ' (res)'] = np.nan
            rowdicts.append(rowdict)
            logging.info('Fit Error = %.1f' % residuals[0, i])

            # leave out the row corresponding with observation 'i'
            logging.info('Cross validation, leaving-one-out: ' + self.obs_ids[i])
            subset = range(n_obs)
            subset.pop(i)
            loo_group_contributions, loo_nullspace = LinearRegression.LeastSquares(
                self.group_matrix[:, subset], self.obs_values[:, subset])
            
            if loo_nullspace.shape[1] > self.group_nullspace.shape[1]:
                logging.warning('example %d is not linearly dependent in the other examples' % i)
                continue
            rowdict['LOO ' + sym + ' (fit)'] = float(np.dot(loo_group_contributions, self.group_matrix[:, i]))
            rowdict['LOO ' + sym + ' (res)'] = \
                rowdict['LOO ' + sym + ' (fit)'] - self.obs_values[0, i]
            logging.info('LOO Error = %.1f' % rowdict['LOO ' + sym + ' (res)'])
        
        logging.info("writing the table of estimation errors for each compound")
        self.html_writer.write('</br><b>Cross-validation table</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.html_writer.write('<font size="1">\n')
        obs_vec = np.matrix([row[sym + ' (obs)'] for row in rowdicts])
        resid_vec = np.matrix([row[sym + ' (res)'] for row in rowdicts])
        rmse = rms_flat(resid_vec.flat)
        
        loo_resid_vec = np.matrix([row['LOO ' + sym + ' (res)']
                                  for row in rowdicts])
        loo_rmse = rms_flat(loo_resid_vec[np.isfinite(loo_resid_vec)].flat)

        self.html_writer.write_ul(['fit RMSE = %.1f [kJ/mol]' % rmse,
                                   'leave-one-out RMSE = %.1f [kJ/mol]' % loo_rmse])
        logging.info("Goodness of fit: RMSE = %.1f [kJ/mol]" % rmse)
        logging.info("Leave-one-out test: RMSE = %.1f [kJ/mol]" % loo_rmse)

        headers = ['Observation',
                   sym + ' (obs)',
                   sym + ' (fit)',
                   sym + ' (res)',
                   'LOO ' + sym + ' (fit)',
                   'LOO ' + sym + ' (res)']
        rowdicts.sort(key=lambda(x):abs(x.get('LOO ' + sym + ' (res.)', 0)),
                      reverse=True)
        self.html_writer.write_table(rowdicts, headers, decimal=1)
        self.html_writer.write('</font>\n')
        self.html_writer.div_end()
        
        self.html_writer.write('</br><b>Cross-validation figure</b>')
        self.html_writer.insert_toggle(start_here=True)
        
        obs_vs_err_fig = plt.figure(figsize=[6.0, 6.0], dpi=100)
        plt.plot(obs_vec.T, resid_vec.T, '.')
        plt.xlabel('Observation')
        plt.ylabel('Estimated (PGC) Residuals')
        plt.hold(True)
        for row in rowdicts:
            if abs(row[sym + ' (res)']) > 2*rmse:
                plt.text(row[sym + ' (obs)'],
                         row[sym + ' (res)'],
                         row['Observation'], fontsize=4,
                         figure=obs_vs_err_fig)
        plt.title('Observed vs. Fitted (PGC) Residuals', figure=obs_vs_err_fig)
        self.html_writer.embed_matplotlib_figure(obs_vs_err_fig)
        self.html_writer.div_end()

    def Mol2Decomposition(self, mol, ignore_protonations=False):
        return self.group_decomposer.Decompose(mol, ignore_protonations, 
                                               strict=True)

    def Reaction2GroupVector(self, sparse):
        total_groupvec = GroupVector(self.groups_data)
        for cid, coeff in sparse.iteritems():
            groupvec = self.cid2groupvec.get(cid, None)
            if groupvec is not None:
                total_groupvec += groupvec * coeff
            else:
                raise MissingCompoundFormationEnergy("C%05d does not have a "
                    "groupvector because - %s" % (cid, self.cid2error[cid]))
        return total_groupvec

    def GetTransfromedReactionEnergies(self, S, cids,
                                       pH=None, I=None, pMg=None, T=None, conc=1):
        pH, I, pMg, T = self.GetConditions(pH=pH, I=I, pMg=pMg, T=T)

        # copy the rows (corresponding to compounds) which are part of the 
        # anchored stoichiometric matrix to a new S_anchored matrix which is
        # in the right order of rows to fit the anchored_contributions vector
        S_anchored = np.matrix(np.zeros((len(self.anchored_cids), S.shape[1])))
        for c, cid in enumerate(cids):
            if cid in self.anchored_cids:
                S_anchored[self.anchored_cids.index(cid), :] = S[c, :]
                
        # calculate the contribution of anchored reaction to the dG0_r of the 
        # desired reactions
        dG0_r = self.anchored_contributions * S_anchored
        
        # normalize out the coefficients of the anchored reactions
        S_anchored = self.anchored_P_L * S_anchored
        S_anchored[np.where(abs(S_anchored) <= self.epsilon)] = 0
        
        # add back the rows of S which correspond to completely new CIDs
        all_cids = list(self.anchored_cids)
        for c, cid in enumerate(cids):
            if cid not in self.anchored_cids:
                all_cids.append(cid)
                S_anchored = np.vstack([S_anchored, S[c, :]])
        
        for r in xrange(S.shape[1]):
            sparse = dict((all_cids[c], S_anchored[c, r])
                          for c in S_anchored[:, r].nonzero()[0].flat)
            try:
                groupvector = self.Reaction2GroupVector(sparse)
                dG0_r[0, r] += self.groupvec2val(groupvector)
            except (GroupMissingTrainDataError, MissingReactionEnergy, MissingCompoundFormationEnergy) as e:
                logging.debug(str(e))
                dG0_r[0, r] = np.nan
        
        ddG0_f = np.matrix(np.zeros((1, S.shape[0])))    
        for c, cid in enumerate(cids):
            diss_table = self.GetDissociationTable(cid)
            if diss_table is not None:
                nH, nMg = self.cid2nH_nMg[cid]
                ddG0_f[0, c] = diss_table.GetDeltaDeltaG0(pH, I, pMg, T, nH=nH, nMg=nMg)
            else:
                ddG0_f[0, c] = np.nan
        
        dG0_r += ddG0_f * S
        
        if conc != 1:
            dG0_r += AddConcentrationsToReactionEnergies(S, cids, T, conc)
        
        return dG0_r

    def get_all_cids(self):
        return sorted(self.cid2groupvec.keys())
    
    def EstimateKeggCids(self):
        """
            Uses the Group Contributions to estimate the entire set of compounds in KEGG,
            and then writes the results to the database as 'gc_pseudoisomers' table
            
            Options:
                override_all_observed_compounds - If True, any observed formation energy is 
                    used instead of the GC estimation. If False, only 'test' compounds are used.
        """
        logging.info("Estimating formation energies for all KEGG")
        
        observed_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='testing')

        for rc in RedoxCarriers().itervalues():
            observed_species.AddPseudoisomer(rc.cid_ox, nH=rc.nH_ox, z=rc.z_ox,
                                             nMg=0, dG0=0.0, ref=rc.ref)
            observed_species.AddPseudoisomer(rc.cid_red, nH=rc.nH_red, z=rc.z_red,
                                             nMg=0, dG0=rc.ddG0, ref=rc.ref)
            observed_species.cid2source_string[rc.cid_ox] = rc.ref
            observed_species.cid2source_string[rc.cid_red] = rc.ref
            
        self.cid2pmap_dict = {}
        self.cid2source_string = {}

        self.html_writer.write('</br><b>Estimated formation energies for KEGG compounds</b>\n')
        self.html_writer.insert_toggle(start_here=True)
        for cid in sorted(self.kegg.get_all_cids()):
            self.html_writer.write('<b>C%05d - %s</b></br>\n' %
                                   (cid, self.kegg.cid2name(cid)))

            diss_table = self.GetDissociationTable(cid)
            if cid in observed_species.get_all_cids():
                pmap_obs = observed_species.cid2PseudoisomerMap(cid)
                self.cid2source_string[cid] = observed_species.cid2SourceString(cid)
                pmatrix = pmap_obs.ToMatrix() # returns a list of (nH, z, nMg, dG0)
                if len(pmatrix) == 1 and diss_table is not None:
                    # assume that only the most abundant pseudoisomer is given
                    # and complete the formation energies of the others using the
                    # pKa values in the dissociation table
                    nH, _z, nMg, dG0 = pmatrix[0]
                    diss_table.SetFormationEnergyByNumHydrogens(dG0=dG0, nH=nH, nMg=nMg)
                    pmap = diss_table.GetPseudoisomerMap()
                    self.SetPseudoisomerMap(cid, pmap)
                else:
                    if diss_table is not None:
                        logging.warning("C%05d has multiple training species, "
                                        "overriding the dissociation table" % cid)
                    self.SetPseudoisomerMap(cid, pmap_obs)
            elif diss_table is None:
                self.html_writer.write('Warning: no dissociation table</br>\n')
                continue
            else:
                nH, nMg = self.cid2nH_nMg[cid]
                groupvector = self.cid2groupvec.get(cid, None)
                if groupvector is None:
                    self.html_writer.write('Warning: no group vector (%s)</br>\n'
                                           % self.cid2error[cid])
                    continue
                try:
                    dG0 = self.groupvec2val(groupvector)
                except GroupMissingTrainDataError as e:
                    # in this case we do not care if a compound violated the group
                    # conservation laws because it might cancel out later when we 
                    # use it to calculate reactions.
                    dG0 = e.value
                    self.html_writer.write('Warning: %s</br>\n' % str(e))
                    logging.debug("C%05d: %s" % (cid, str(e)))

                self.cid2source_string[cid] = "Group Contribution"
                
                if self.transformed:
                    diss_table.SetTransformedFormationEnergy(dG0_tag=dG0, 
                        pH=default_pH, I=default_I, pMg=default_pMg, T=default_T)
                else:
                    diss_table.SetFormationEnergyByNumHydrogens(dG0=dG0,
                                                                nH=nH, nMg=nMg)
                pmap = diss_table.GetPseudoisomerMap()
                self.SetPseudoisomerMap(cid, pmap)
            self.cid2PseudoisomerMap(cid).WriteToHTML(self.html_writer)

        self.html_writer.div_end()
        
        logging.info("Writing the results to the database")
        self.ToDatabase(self.db, self.THERMODYNAMICS_TABLE_NAME)

    def VerifyReaction(self, sparse):
        """
            Inherited from Thermodynamics.
            In this case, should check that the total reaction groupvec is orthogonal to 
            the kernel.
        """
        try:
            groupvec = self.Reaction2GroupVector(sparse)
        except MissingCompoundFormationEnergy as e:
            raise MissingReactionEnergy(str(e), sparse)

        try:
            # this will practically check that the overall
            # reaction groupvec is orthogonal to the kernel
            self.groupvec2val(groupvec)
        except GroupMissingTrainDataError as e:
            raise MissingReactionEnergy(e.Explain(self), sparse)
        
    def EstimateKeggRids(self, pH=None, pMg=None, I=None ,T=None):
        counters = {'ok':0, 'balance':0, 'formation':0, 'reaction':0, 'key':0}
        for rid in sorted(self.kegg.get_all_rids()):
            reaction = self.kegg.rid2reaction(rid)
            try:
                reaction.Balance(balance_water=True)
                dG0_r = reaction.PredictReactionEnergy(self)
                print "R%05d = %7.2f" % (rid, dG0_r)
                counters['ok'] += 1
            except KeggReactionNotBalancedException as e:
                print "R%05d: KeggReactionNotBalancedException %s" % (rid, str(e))
                counters['balance'] += 1
            except MissingCompoundFormationEnergy as e:
                print "R%05d: Missing the formation energy of C%05d (%s)" % (rid, e.cid, self.kegg.cid2name(e.cid))
                counters['formation'] += 1
            except MissingReactionEnergy as e:
                print "R%05d: MissingReactionEnergy %s" % (rid, str(e))
                counters['reaction'] += 1
            except KeyError as e:
                print "R%05d: KeyError %s" % (rid, str(e))
                counters['key'] += 1
        
        for k, v in counters.iteritems():
            print k, v
            
    def cid2PseudoisomerMap(self, cid):
        """
            Overrides the cid2PseudoisomerMap method from PsuedoisomerTableThermodynamics
            in order to raise the group-contribution-specific error using the 
            MissingCompoundFormationEnergy exception.
        """
        if not self.cid2pmap_dict:
            raise Exception("You must run the method EstimateKeggCids before using this one.")
        elif cid in self.cid2pmap_dict and self.cid2pmap_dict[cid] is not None:
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy(self.cid2error.get(cid, ""), cid)
        
    def GroupMatrixRowToString(self, r):
        nonzero_columns = list(np.where(abs(r) > self.epsilon)[0].flat)
        return " | ".join(["%g : %s" % (r[j], self.groups_data.all_groups[j].name)
                           for j in nonzero_columns.flat])
            
    def GetGroupContribution(self, name, nH, z, nMg=0):
        gr = Group(None, name, nH, z, nMg)
        gid = self.groups_data.Index(gr)
        return self.group_contributions[0, gid]
    
    def GetPkaOfGroup(self, name, nH, z, nMg=0):
        dG0_gr_deprotonated = self.GetGroupContribution(name, nH-1, z-1, nMg)
        dG0_gr_protonated = self.GetGroupContribution(name, nH, z, nMg)
        if not dG0_gr_deprotonated or not dG0_gr_protonated:
            return None
        pKa = (dG0_gr_deprotonated - dG0_gr_protonated)/(R*default_T*np.log(10))
        return pKa
    
    def GetPkmgOfGroup(self, name, nH, z, nMg):
        dG0_gr_without_Mg = self.GetGroupContribution(name, nH, z-2, nMg-1)
        dG0_gr_with_Mg = self.GetGroupContribution(name, nH, z, nMg)
        if not dG0_gr_without_Mg or not dG0_gr_with_Mg:
            return None
        pK_Mg = (dG0_gr_without_Mg + dG0_f_Mg - dG0_gr_with_Mg)/(R*default_T*np.log(10))
        return pK_Mg

    def AnalyzeSingleKeggCompound(self, cid, ignore_protonations=False):
        print 'Analyzing C%05d (%s):' % (cid, self.kegg.cid2name(cid))
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            print "This compounds doesn't have a dissociation table"
            return
        mol = diss_table.GetMostAbundantMol(pH=default_pH, I=0, 
                                            pMg=14, T=default_T)
        if mol is None:
            print "This compounds dissociation table does not have a Mol description"
            return

        try:
            decomposition = self.Mol2Decomposition(mol,
                                                   ignore_protonations=ignore_protonations)
            print decomposition.ToTableString()
            print 'nH =', decomposition.Hydrogens()
            print 'z =', decomposition.NetCharge()
            print 'nMg = ', decomposition.Magnesiums()
        except GroupDecompositionError as e:
            print "Cannot decompose compound to groups: " + str(e)
        mol.Draw()
        
    def AnalyzeSingleCompound(self, mol):
        decomposition = self.group_decomposer.Decompose(mol, 
                                ignore_protonations=False, strict=False)
        print decomposition.ToTableString()
        print 'nH =', decomposition.Hydrogens()
        print 'z =', decomposition.NetCharge()
        print 'nMg = ', decomposition.Magnesiums()
        mol.Draw()
        
    def init(self):
        self.LoadGroups(True)
        self.LoadObservations(True)
        self.LoadGroupVectors(True)
        if self.db.DoesTableExist(self.CONTRIBUTION_TABLE_NAME):
            self.LoadContributionsFromDB()
        else:
            self.Train()
            self.EstimateKeggCids()
        
        reader = self.db.DictReader(self.THERMODYNAMICS_TABLE_NAME)
        PsuedoisomerTableThermodynamics._FromDictReader(
            reader, self, label=None, name="Group Contribution",
            warn_for_conflicting_refs=False)

#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-b", "--biochemical", action="store_true",
                          dest="transformed", default=False,
                          help="Use biochemical (transformed) Group Contributions")
    opt_parser.add_option("-t", "--train", action="store_true",
                          dest="train_only", default=False,
                          help="A flag for running the TRAIN only (without TEST)")
    opt_parser.add_option("-e", "--test", action="store_true",
                          dest="test_only", default=False,
                          help="A flag for running the TEST only (without TRAIN)")
    opt_parser.add_option("-d", "--from_database", action="store_true",
                          dest="from_database", default=False,
                          help="A flag for loading the data from the DB instead of "
                               "the CSV files (saves time but no debug information)")
    return opt_parser

if __name__ == '__main__':
    options, _ = MakeOpts().parse_args(sys.argv)
    util._mkdir('../res')
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    
    # use the flag -i or --train for train only
    # use the flag -e or --test for test only
    if options.transformed:
        prefix = 'bgc'
    else:
        prefix = 'pgc'
    
    if options.test_only:
        html_writer = HtmlWriter('../res/%s_test.html' % prefix)
    elif options.train_only:
        html_writer = HtmlWriter('../res/%s_train.html' % prefix)
    else:
        html_writer = HtmlWriter('../res/%s.html' % prefix)
        
    G = GroupContribution(db=db, html_writer=html_writer,
                          transformed=options.transformed)
    
    G.LoadGroups(options.from_database)
    G.LoadObservations(options.from_database)
    G.LoadGroupVectors(options.from_database)
    
    if options.test_only:
        G.LoadContributionsFromDB()
    else:
        G.Train()
        G.WriteRegressionReport()
        G.AnalyzeTrainingSet()
    
    if not options.train_only:
        G.EstimateKeggCids()

