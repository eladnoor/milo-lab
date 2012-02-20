#!/usr/bin/python

import logging, types, json, sys
from collections import defaultdict

from copy import deepcopy
import pylab
import numpy as np
from pygibbs.thermodynamic_constants import R, default_pH, default_T,\
    dG0_f_Mg, default_I, default_pMg, RedoxCarriers
from pygibbs.thermodynamics import MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.groups_data import Group, GroupsData
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.molecule import Molecule
from pygibbs.group_observation import GroupObervationCollection
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
        for i in self.kernel_rows:
            nonzero_columns = np.where(abs(gc.group_nullspace[i, :]) > 1e-10)[0]
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
        if not FromDatabase:
            fname = "../data/thermodynamics/groups_species.csv"
            self.groups_data = GroupsData.FromGroupsFile(fname,
                                                         transformed=self.transformed)
            self.groups_data.ToDatabase(self.db)
            self.group_decomposer = GroupDecomposer(self.groups_data)
        else:
            self.groups_data = GroupsData.FromDatabase(self.db,
                                                       transformed=self.transformed)
            self.group_decomposer = GroupDecomposer(self.groups_data)

    def LoadObservations(self, FromDatabase=False):
        if not FromDatabase:
            logging.info("Reading observations from files")
            dissociation = self.GetDissociationConstants()
            self.obs_collection = GroupObervationCollection.FromFiles(
                                    html_writer=self.html_writer, 
                                    dissociation=dissociation,
                                    transformed=self.transformed)
            self.obs_collection.ToDatabase(self.db, self.OBSERVATION_TABLE_NAME)
            self.cid2nH_nMg = self.obs_collection.cid2nH_nMg
        else:
            logging.info("Reading observations from database")
            self.obs_collection = GroupObervationCollection.FromDatabase(
                                    db=self.db,
                                    table_name=self.OBSERVATION_TABLE_NAME,
                                    transformed=self.transformed)
        
        self.obs_collection.ReportToHTML()

    def LoadGroupVectors(self, FromDatabase=False):
        self.cid2groupvec = {}
        self.cid2error = {}            

        if not FromDatabase:
            logging.info("Decomposing all compounds and calculating group vectors")
            self.html_writer.write('</br><b>All Groupvectors</b>\n')
            self.html_writer.insert_toggle(start_here=True)
            # When using non-transformed energies, it is very important for
            # the group vectors of each compound to represent the correct
            # pseudoisomer (same nH and nMg used in the reverse transform and
            # the list of formation energies). Here we use the dictionary 
            # self.cid2nH_nMg that is copied from GroupObervationCollection
            
            dissociation = self.GetDissociationConstants()
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
        else:
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

    def Train(self):
        logging.info("Calculating the linear regression data")
        cids, S, b, anchored = self.obs_collection.GetStoichiometry()
        obs_ids = [obs.obs_id for obs in self.obs_collection.observations]
        obs_types = [obs.obs_type for obs in self.obs_collection.observations]

        anchored_cols = np.where(anchored==1)[1]
        # now remove anchored data from S and leave only the data which will be 
        # used for calculating the group contributions
        g, _ = LinearRegression.LeastSquares(S[:, anchored_cols],
                                             b[:, anchored_cols])
        P_C, P_L = LinearRegression.ColumnProjection(S[:, anchored_cols])
        b -= np.dot(np.dot(g, P_C), S)
        S = np.dot(P_L, S)
        
        # set epsilon-small values to absolute 0
        S[np.where(abs(S) < 1e-10)] = 0
        
        # removed zero rows (compounds) from S
        used_cid_indices = set(np.nonzero(np.sum(abs(S), 1))[0])
        for i_cid, cid in enumerate(cids):
            if self.cid2groupvec[cid] is None:
                used_cid_indices.difference_update([i_cid])
                for i_obs in np.nonzero(S[i_cid, :])[0]:
                    logging.warning("%s is removed because C%05d has no group vector, "
                                    "but is still part of the final stoichiometric matrix"
                                    % (obs_ids[i_obs], cid))
                    S[:, i_obs] = 0

        used_cid_indices = sorted(used_cid_indices)
        S = S[used_cid_indices, :]

        # removed zero column (observations) from S
        #nonzero_cols = np.nonzero(np.sum(abs(S), 0))[0]
        #nonzero_cols = set(nonzero_cols).difference(set_of_obs_to_remove)
        #nonzero_cols = sorted(nonzero_cols)
        #S = S[:, nonzero_cols]
        #b = b[:, nonzero_cols]
        #observations = [observations[i] for i in nonzero_cols]
        
        n_groups = len(self.groups_data.GetGroupNames()) # number of groups
        G = np.zeros((len(used_cid_indices), n_groups))
        for i, i_cid in enumerate(used_cid_indices):
            G[i, :] = self.cid2groupvec[cids[i_cid]].Flatten()

        GS = np.dot(G.T, S)

        # 'unique' the rows GS. For each set of rows that is united,
        # the Y-value for the new row is the average of the corresponding Y-values.
        unique_GS, col_mapping = LinearRegression.ColumnUnique(GS, remove_zero=True)
        unique_b = np.zeros((1, unique_GS.shape[1]))
        unique_obs_types = []
        unique_obs_ids = []
        for i, old_indices in sorted(col_mapping.iteritems()):
            unique_b[0, i] = np.mean(b[0, old_indices])
            unique_obs_types.append(obs_types[old_indices[0]]) # take the type of the first one (not perfect...)
            unique_obs_ids.append(', '.join([obs_ids[i] for i in old_indices]))            
        
        self.group_matrix = unique_GS
        self.obs_values = unique_b
        self.obs_ids = unique_obs_ids
        self.obs_types = unique_obs_types

        logging.info("Performing linear regression")
        self.group_contributions, self.group_nullspace = \
            LinearRegression.LeastSquares(self.group_matrix, self.obs_values,
                                          reduced_row_echlon=False)
        
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
            
    def does_table_exist(self, table_name):
        for unused_ in self.db.Execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def groupvec2val(self, groupvec):
        if self.group_contributions == None or self.group_nullspace == None:
            raise Exception("You need to first Train the system before using it to estimate values")

        gv = np.array(groupvec.Flatten()).T
        val = float(np.dot(self.group_contributions, gv))
        v = abs(np.dot(self.group_nullspace.T, gv))
        k_list = [i for i in np.where(v > 1e-10)[0]]
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
        for j, dG0_gr in enumerate(self.group_contributions[0, :]):
            obs_lists_dict = defaultdict(list)
            for k in self.group_matrix[j, :].nonzero()[0]:
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
        deviations = [] # a list of all the data about every example
        loo_deviations = [] # a list of only the examples which are included in the LOO test
        
        for i in xrange(n_obs):
            row = {'name':self.obs_ids[i], 'obs':self.obs_values[0, i]}
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.obs_types[i] not in ['formation', 'reaction']:
                continue
            
            row['fit'] = float(np.dot(self.group_contributions, self.group_matrix[:, i]))
            row['fit_resid'] = row['fit']-row['obs'] 
            deviations.append(row)
            logging.info('Fit Error = %.1f' % (row['fit_resid']))

            # leave out the row corresponding with observation 'i'
            logging.info('Cross validation, leaving-one-out: ' + self.obs_ids[i])
            subset = range(n_obs)
            subset.pop(i)
            loo_group_contributions, loo_nullspace = LinearRegression.LeastSquares(
                self.group_matrix[:, subset], self.obs_values[:, subset])
            
            if loo_nullspace.shape[1] > self.group_nullspace.shape[1]:
                logging.warning('example %d is not linearly dependent in the other examples' % i)
                continue
            row['loo'] = float(np.dot(loo_group_contributions, self.group_matrix[:, i]))
            row['loo_resid'] = row['loo'] - row['obs']
            loo_deviations.append(row)
            logging.info('LOO Error = %.1f' % row['loo_resid'])
        
        logging.info("writing the table of estimation errors for each compound")
        self.html_writer.write('</br><b>Cross-validation table</b>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)

        obs_vec = np.array([row['obs'] for row in deviations])
        fit_vec = np.array([row['fit'] for row in deviations])
        
        resid_vec = np.array([row['fit_resid'] for row in deviations])
        rmse = pylab.rms_flat(resid_vec)
        
        loo_resid_vec = np.array([row['loo_resid'] for row in loo_deviations])
        loo_rmse = pylab.rms_flat(loo_resid_vec)

        self.html_writer.write_ul(['fit_rmse(pred) = %.2f kJ/mol' % rmse,
                                   'loo_rmse(pred) = %.2f kJ/mol' % loo_rmse])
        logging.info("Goodness of fit: RMSE = %.2f kJ/mol" % rmse)
        logging.info("Leave-one-out test: RMSE = %.2f kJ/mol" % loo_rmse)

        self.html_writer.write('<font size="1">\n')
        self.html_writer.table_start()
        headers = ['Observation Name',
                   '&Delta;<sub>f</sub>G<sub>obs</sub> [kJ/mol]',
                   '&Delta;<sub>f</sub>G<sub>fit</sub> [kJ/mol]',
                   'Residual [kJ/mol]',
                   '&Delta;<sub>f</sub>G<sub>LOO</sub> [kJ/mol]',
                   'Leave-one-out Residual [kJ/mol]']
        self.html_writer.table_writerow(headers)
        deviations.sort(key=lambda(x):abs(x.get('loo_resid', 0)), reverse=True)
        for row in deviations:
            table_row = [row['name']]
            table_row += ['%.1f' % row[key] for key in ['obs', 'fit', 'fit_resid']]
            if 'loo' in row:
                table_row += ['%.1f' % row[key] for key in ['loo', 'loo_resid']]
            else:
                table_row += ['N/A', 'N/A']
            self.html_writer.table_writerow(table_row)
        self.html_writer.table_end()
        self.html_writer.write('</font>')
        self.html_writer.div_end()
        
        logging.info("Plotting graphs for PGC vs. observed data")
        self.html_writer.write('</br><b>Cross-validation figure 1</b>')
        self.html_writer.insert_toggle(start_here=True)

        obs_vs_est_fig = pylab.figure(figsize=[6.0, 6.0], dpi=100)
        pylab.plot(obs_vec, fit_vec, '.', figure=obs_vs_est_fig)
        pylab.xlabel('Observation', figure=obs_vs_est_fig)
        pylab.ylabel('Estimation (PGC)', figure=obs_vs_est_fig)
        pylab.hold(True)
        for row in deviations:
            if abs(row['fit_resid']) > 2*rmse:
                pylab.text(row['obs'], row['fit'], row['name'], fontsize=4,
                           figure=obs_vs_est_fig)
        pylab.title('Observed vs. Estimated (PGC)', figure=obs_vs_est_fig)
        self.html_writer.embed_matplotlib_figure(obs_vs_est_fig)
        self.html_writer.div_end()

        self.html_writer.write('</br><b>Cross-validation figure 2</b>')
        self.html_writer.insert_toggle(start_here=True)
        
        obs_vs_err_fig = pylab.figure(figsize=[6.0, 6.0], dpi=100)
        pylab.plot(obs_vec, resid_vec, '.')
        pylab.xlabel('Observation')
        pylab.ylabel('Estimated (PGC) Residuals')
        pylab.hold(True)
        for row in deviations:
            if abs(row['fit_resid']) > 2*rmse:
                pylab.text(row['obs'], row['fit_resid'], row['name'], fontsize=4)
        pylab.title('Observed vs. Estimated (PGC) Residuals', figure=obs_vs_est_fig)
        self.html_writer.embed_matplotlib_figure(obs_vs_err_fig)
        self.html_writer.div_end()

    def Mol2Decomposition(self, mol, ignore_protonations=False):
        return self.group_decomposer.Decompose(mol, ignore_protonations, 
                                               strict=True)

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

        self.html_writer.write('</br><b>Estimated formation energies for KEGG compounds</b></br>\n')
        self.html_writer.insert_toggle(start_here=True)
        for cid in sorted(self.kegg.get_all_cids()):
            self.html_writer.write('<b>C%05d - %s</b></br>\n' %
                                   (cid, self.kegg.cid2name(cid)))

            diss_table = self.GetDissociationTable(cid)
            pmap = None
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

    def Reaction2GroupVector(self, sparse):
        total_groupvec = GroupVector(self.groups_data)
        for cid, coeff in sparse.iteritems():
            if cid not in self.cid2groupvec:
                # one of the compounds has no observed or estimated formation energy
                raise MissingReactionEnergy("%s (C%05d) has no GroupVector" % 
                    (self.kegg.cid2name(cid), cid), sparse)
            groupvec = self.cid2groupvec[cid]
            if groupvec is not None:
                total_groupvec += groupvec * coeff
        return total_groupvec
    
    def VerifyReaction(self, sparse):
        """
            Inherited from Thermodynamics.
            In this case, should check that the total reaction groupvec is orthogonal to 
            the kernel.
        """
        groupvec = self.Reaction2GroupVector(sparse)
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
        nonzero_columns = np.where(abs(r) > 1e-10)[0]
        return " | ".join(["%g : %s" % (r[j], self.groups_data.all_groups[j].name)
                           for j in nonzero_columns])
            
    def LoadContributionsFromDB(self):
        logging.info("loading the group contribution data from the database")
        self.group_contributions = []
        for row in self.db.DictReader(self.CONTRIBUTION_TABLE_NAME):
            self.group_contributions.append(row['dG0_gr'])
        self.group_contributions = pylab.array([self.group_contributions])
        self.group_nullspace = self.db.LoadSparseNumpyMatrix(self.NULLSPACE_TABLE_NAME)
                        
    def GetGroupContribution(self, name, nH, z, nMg=0):
        gr = Group(None, name, nH, z, nMg)
        gid = self.groups_data.Index(gr)
        return self.group_contributions[0, gid]
    
    def GetPkaOfGroup(self, name, nH, z, nMg=0):
        dG0_gr_deprotonated = self.GetGroupContribution(name, nH-1, z-1, nMg)
        dG0_gr_protonated = self.GetGroupContribution(name, nH, z, nMg)
        if not dG0_gr_deprotonated or not dG0_gr_protonated:
            return None
        pKa = (dG0_gr_deprotonated - dG0_gr_protonated)/(R*default_T*pylab.log(10))
        return pKa
    
    def GetPkmgOfGroup(self, name, nH, z, nMg):
        dG0_gr_without_Mg = self.GetGroupContribution(name, nH, z-2, nMg-1)
        dG0_gr_with_Mg = self.GetGroupContribution(name, nH, z, nMg)
        if not dG0_gr_without_Mg or not dG0_gr_with_Mg:
            return None
        pK_Mg = (dG0_gr_without_Mg + dG0_f_Mg - dG0_gr_with_Mg)/(R*default_T*pylab.log(10))
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
        self.LoadContributionsFromDB()
        
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
    opt_parser.add_option("-c", "--compound", action="store", type="int",
                          dest="cid", default=None,
                          help="The KEGG ID of a compound")
    opt_parser.add_option("-r", "--reaction", action="store", type="int",
                          dest="rid", default=None,
                          help="The KEGG ID of a reaction")
    opt_parser.add_option("-s", "--smiles", action="store", type="string",
                          dest="smiles", default=None,
                          help="A SMILES string of a compound")
    opt_parser.add_option("-i", "--inchi", action="store", type="string",
                          dest="inchi", default=None,
                          help="An InChI string of a compound")
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
    
    if options.smiles or options.inchi or options.cid or options.rid:
        G = GroupContribution(db=db)
        if options.smiles: # -s <SMILES>
            print 'Analyzing SMILES %s:' % (options.smiles)
            mol = Molecule._FromFormat(options.smiles, 'smiles')
            G.LoadGroupsFromFile()
            G.AnalyzeSingleCompound(mol)
        elif options.inchi: #-i <INCHI>
            print 'Analyzing InChI %s:' % (options.inchi)
            mol = Molecule._FromFormat(options.inchi, 'inchi')
            G.LoadGroupsFromFile()
            G.AnalyzeSingleCompound(mol)
        elif options.cid: # -c <CID>
            print 'Analyzing Compound C%05d:' % (options.cid)
            G.init()
            G.AnalyzeSingleKeggCompound(options.cid, ignore_protonations=True)
        elif options.rid: # -r <RID>
            print 'Analyzing Reaction R%05d:' % (options.rid)
            G.init()
            reaction = G.kegg.rid2reaction(options.rid)
            dG0_r = reaction.PredictReactionEnergy(G)
            print "dG0_r = %.2f" % dG0_r
    else:
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
        
        if not options.test_only:
            G.LoadGroups(options.from_database)
            G.LoadObservations(options.from_database)
            G.LoadGroupVectors(options.from_database)
            G.Train()
            G.WriteRegressionReport()
            G.AnalyzeTrainingSet()
        else:
            G.LoadGroups(True)
            G.LoadObservations(True)
            G.LoadGroupVectors(True)
            G.LoadContributionsFromDB()
        
        if not options.train_only:
            G.EstimateKeggCids()

