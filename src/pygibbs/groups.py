#!/usr/bin/python

import csv, logging, types, json, sys

from copy import deepcopy
import pylab
import numpy as np
from pygibbs.thermodynamic_constants import R, default_pH, default_T,\
    dG0_f_Mg, default_I, default_pMg, F, RedoxCarriers
from pygibbs.thermodynamics import MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.groups_data import Group, GroupsData
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.sparse_kernel import SparseKernel, CplexNotInstalledError
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
        self.sparse_kernel = False

        self.group_nullspace = None
        self.group_contributions = None
        self.obs_collection = None
        
        self.cid2error = {}
        self.cid2groupvec = None

        if transformed:
            prefix = 'bgc'
        else:
            prefix = 'pgc'
        
        self.PMAP_TABLE_NAME = prefix + '_pseudoisomers'
        self.ERROR_TABLE_NAME = prefix + '_errors'
        self.GROUPVEC_TABLE_NAME = prefix + '_groupvector'
        self.NULLSPACE_TABLE_NAME = prefix + '_nullspace'
        self.CONTRIBUTION_TABLE_NAME = prefix + '_contribution'
        self.REGRESSION_TABLE_NAME = prefix + '_regression'
        
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
    
    def init(self):
        self.LoadGroupsFromDatabase()
        self.LoadContributionsFromDB()
        if self.db.DoesTableExist(self.PMAP_TABLE_NAME):
            self.FromDatabase(self.db, table_name=self.PMAP_TABLE_NAME)
        
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
        
    def LoadGroupsFromFile(self):
        #if self.transformed:
        #    fname = "../data/thermodynamics/groups_species_transformed.csv"
        #else:
        fname = "../data/thermodynamics/groups_species.csv"
        self.groups_data = GroupsData.FromGroupsFile(fname,
                                                     transformed=self.transformed)
        self.groups_data.ToDatabase(self.db)
        self.group_decomposer = GroupDecomposer(self.groups_data)
    
    def LoadGroupsFromDatabase(self):
        self.groups_data = GroupsData.FromDatabase(self.db,
                                                   transformed=self.transformed)
        self.group_decomposer = GroupDecomposer(self.groups_data)
            
    def Train(self, FromFiles=True):
        self.obs_collection = GroupObervationCollection(
                                db=self.db, 
                                html_writer=self.html_writer, 
                                group_decomposer=self.group_decomposer, 
                                dissociation=self.GetDissociationConstants(),
                                transformed=self.transformed)
        if FromFiles:
            self.obs_collection.FromFiles()
            self.obs_collection.ToDatabase()
            self.obs_collection.ToCSV('../res/gc_pseudoisomers.csv',
                                      '../res/gc_observations.csv')
        else:
            self.obs_collection.FromDatabase()
        
        self.group_matrix, self.obs_values, self.obs_types, self.obs_names = \
            self.obs_collection.GetRegressionData()
        self.SaveRegressionDataToDB()

        self.group_contributions, self.group_nullspace = self.RunLinearRegression()
        self.SaveContributionsToDB()

    def SaveRegressionDataToDB(self):
        self.db.CreateTable(self.REGRESSION_TABLE_NAME,
                            'name TEXT, type TEXT, groupvec TEXT, dG0 REAL')
        for r in xrange(self.group_matrix.shape[0]):
            groupvec = GroupVector(self.groups_data, self.group_matrix[r, :])
            self.db.Insert(self.REGRESSION_TABLE_NAME,
                           [self.obs_names[r], self.obs_types[r], 
                            groupvec.ToJSONString(), self.obs_values[r, 0]])

    def SaveContributionsToDB(self):
        logging.info("storing the group contribution data in the database")
        
        # write a table of the group contributions
        self.db.CreateTable(self.CONTRIBUTION_TABLE_NAME,
                            'gid INT, name TEXT, protons INT, charge INT, '
                            'nMg INT, dG0_gr REAL, nullspace TEXT')
        for j, group_name in enumerate(self.groups_data.GetGroupNames()):
            dG0_gr = self.group_contributions[j, 0]
            nullspace_str = ','.join(["%.2f" % x for x in self.group_nullspace[:, j]])
            
            if self.transformed:
                nH, z, nMg = None, None, None
            else:
                group = self.groups_data.all_groups[j]
                nH, z, nMg = group.hydrogens, group.charge, group.nMg

            self.db.Insert(self.CONTRIBUTION_TABLE_NAME,
                [j, group_name, nH, z, nMg, dG0_gr, nullspace_str])
            
        self.db.CreateTable(self.NULLSPACE_TABLE_NAME, 'dimension INT, group_vector TEXT')
        for i in xrange(self.group_nullspace.shape[0]):
            groupvec = GroupVector(self.groups_data, self.group_nullspace[i, :])
            groupvec.RemoveEpsilonValues(epsilon=1e-10)
            self.db.Insert(self.NULLSPACE_TABLE_NAME, 
                           [i, groupvec.ToJSONString()])

        self.db.Commit()
            
    def RunLinearRegression(self):
        group_contributions, nullspace = LinearRegression.LeastSquares(
                        self.group_matrix, self.obs_values, reduced_row_echlon=False)
        
        if self.sparse_kernel:
            try:
                nullspace = SparseKernel(self.group_matrix).Solve()
            except CplexNotInstalledError:
                logging.warning("CPLEX is not installed on this system, using a non-sparse"
                                " method for describing the Kernel of the group matrix")
            
        return group_contributions, nullspace
    
    def does_table_exist(self, table_name):
        for unused_ in self.db.Execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def groupvec2val(self, groupvec):
        if self.group_contributions == None or self.group_nullspace == None:
            raise Exception("You need to first Train the system before using it to estimate values")

        gv = np.array(groupvec.Flatten())
        val = float(np.dot(gv, self.group_contributions))
        v = abs(np.dot(self.group_nullspace, gv))
        k_list = [i for i in np.where(v > 1e-10)[0]]
        if k_list:
            raise GroupMissingTrainDataError(val, "can't estimate because the input "
                "is not orthogonal to the kernel", k_list)
        return val
    
    def WriteRegressionReport(self, T=default_T, pH=default_pH):        
        self.html_writer.write('</br><b>Regression report</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.html_writer.write_ul(['observations: %d' % self.group_matrix.shape[0],
           'groups: %d' % self.group_matrix.shape[1],
           'rank: %d' % LinearRegression.MatrixRank(self.group_matrix)])
        self.html_writer.write('</br><table border="1">\n<tr>'
            '<td width="5%%">#</td>'
            '<td width="20%%">Name</td>'
            '<td width="5%%">&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td>'
            '<td width="70%%">Group Vector</td></tr>')
        for i in xrange(self.group_matrix.shape[0]):
            groupvec = GroupVector(self.groups_data, self.group_matrix[i, :])
            self.html_writer.write('<tr>' +
                '<td width="5%%">%d</td>' % i +
                '<td width="20%%">%s</td>' % self.obs_names[i] + 
                '<td width="5%%">%8.2f</td>' % self.obs_values[i] +
                '<td width="70%%">%s</td></tr>' % str(groupvec))
        self.html_writer.write('</table>')
        self.html_writer.div_end()
        
        self.html_writer.write('</br><b>Group Contributions</b>\n')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)
        self.html_writer.write('</br><font size="1">\n')
        dict_list = []
        group_names = self.groups_data.GetGroupNames()
        for j, dG0_gr in enumerate(self.group_contributions):
            obs_lists_dict = {'acid-base':[], 'formation':[], 'reaction':[]}
            for k in self.group_matrix[:, j].nonzero()[0]:
                obs_lists_dict[self.obs_types[k]].append(self.obs_names[k])
            d = {"#":"%d" % j, "Group Name":group_names[j], 
                 "&#x394;<sub>gr</sub>G [kJ/mol]":"%8.2f" % dG0_gr,
                 "dissociations":' | '.join(obs_lists_dict['acid-base']),
                 "formations":' | '.join(obs_lists_dict['formation']),
                 "reactions":' | '.join(obs_lists_dict['reaction'])}
            if not self.transformed:
                group = self.groups_data.all_groups[j]
                d["nH"] = group.hydrogens
                d["charge"] = group.charge
                d["nMg"] = group.nMg
            dict_list.append(d)
        if not self.transformed:
            self.html_writer.write_table(dict_list, headers=["#", "Group Name",
                "nH", "charge", "nMg", "&#x394;<sub>gr</sub>G [kJ/mol]", 
                "dissociations", "formations", "reactions"])
        else:
            self.html_writer.write_table(dict_list, headers=["#", "Group Name",
                "&#x394;<sub>gr</sub>G [kJ/mol]", 
                "dissociations", "formations", "reactions"])
        self.html_writer.write('</font>\n')
        self.html_writer.div_end()

        # Null-space matrix
        self.html_writer.write('</br><b>Nullspace of regression matrix</b>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)
        self.html_writer.write('</br>')
        self.db.Table2HTML(self.html_writer, self.NULLSPACE_TABLE_NAME)
        self.html_writer.div_end()

    def analyze_training_set(self):
        n_obs = self.group_matrix.shape[0]
        deviations = [] # a list of all the data about every example
        loo_deviations = [] # a list of only the examples which are included in the LOO test
        
        for i in xrange(n_obs):
            row = {'name':self.obs_names[i], 'obs':self.obs_values[i, 0]}
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.obs_types[i] not in ['formation', 'reaction']:
                continue
            
            row['fit'] = float(np.dot(self.group_matrix[i, :], self.group_contributions))
            row['fit_resid'] = row['fit']-row['obs'] 
            deviations.append(row)
            logging.info('Fit Error = %.1f' % (row['fit_resid']))

            # leave out the row corresponding with observation 'i'
            logging.info('Cross validation, leaving-one-out: ' + self.obs_names[i])
            subset = range(n_obs)
            subset.pop(i)
            loo_group_contributions, loo_nullspace = LinearRegression.LeastSquares(
                self.group_matrix[subset, :], self.obs_values[subset, :])
            
            if loo_nullspace.shape[0] > self.group_nullspace.shape[0]:
                logging.warning('example %d is not linearly dependent in the other examples' % i)
                continue
            row['loo'] = float(np.dot(self.group_matrix[i, :], loo_group_contributions))
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
                   '&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]',
                   '&#x394;<sub>f</sub>G<sub>fit</sub> [kJ/mol]',
                   'Residual [kJ/mol]',
                   '&#x394;<sub>f</sub>G<sub>LOO</sub> [kJ/mol]',
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
        all_cids = set(self.cid2pmap_dict.keys() + self.cid2error.keys())
        return sorted(all_cids)
    
    def EstimateKeggCids(self, cid_list=None):
        """
            Uses the Group Contributions to estimate the entire set of compounds in KEGG,
            and then writes the results to the database as 'gc_pseudoisomers' table
            
            Options:
                override_all_observed_compounds - If True, any observed formation energy is 
                    used instead of the GC estimation. If False, only 'test' compounds are used.
        """
        logging.info("Estimating formation energies for all KEGG. Please be patient for a few minutes...")
        
        observed_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='testing')

        for rc in RedoxCarriers().itervalues():
            observed_species.AddPseudoisomer(rc.cid_ox, nH=rc.nH_ox, z=rc.z_ox,
                                             nMg=0, dG0=0.0, ref=rc.ref)
            observed_species.AddPseudoisomer(rc.cid_red, nH=rc.nH_red, z=rc.z_red,
                                             nMg=0, dG0=rc.ddG0, ref=rc.ref)
            observed_species.cid2SourceString[rc.cid_ox] = rc.ref
            observed_species.cid2SourceString[rc.cid_red] = rc.ref
            
        self.cid2pmap_dict = {}
        self.cid2source_string = {}
        self.cid2error = {}
        self.cid2groupvec = {}

        self.html_writer.write('</br><b>Estimated formation energies for KEGG compounds</b></br>\n')
        self.html_writer.insert_toggle(start_here=True)
        for cid in sorted(cid_list or self.kegg.get_all_cids()):
            self.html_writer.write('<b>C%05d - %s</b></br>\n' % (cid, self.kegg.cid2name(cid)))

            diss_table = self.GetDissociationTable(cid)
            pmap = None
            if cid in observed_species.get_all_cids():
                pmap_obs = observed_species.cid2PseudoisomerMap(cid)
                groupvector = GroupVector(self.groups_data) # use an empty group vector
                source_string = observed_species.cid2SourceString(cid)

                pmatrix = pmap_obs.ToMatrix() # returns a list of (nH, z, nMg, dG0)
                if len(pmatrix) == 1 and diss_table is not None:
                    # assume that only the most abundant pseudoisomer is given
                    # and complete the formation energies of the others using the
                    # pKa values in the dissociation table
                    nH, _z, nMg, dG0 = pmatrix[0]
                    diss_table.SetFormationEnergyByNumHydrogens(dG0=dG0, nH=nH, nMg=nMg)
                    pmap = diss_table.GetPseudoisomerMap()
                else:
                    if diss_table is not None:
                        logging.warning("C%05d has multiple training species, "
                                        "overriding the dissociation table" % cid)
                    pmap = pmap_obs
            elif diss_table is None:
                self.html_writer.write('Warning: no dissociation table</br>\n')
                continue
            else:
                nH, nMg = diss_table.GetMostAbundantPseudoisomer(
                    pH=default_pH, I=0, pMg=14, T=default_T)
                mol = diss_table.GetMol(nH, nMg)
                if mol is None:
                    continue
                try:
                    decomposition = self.Mol2Decomposition(mol,
                                                           ignore_protonations=False)
                except GroupDecompositionError as e:
                    self.cid2error[cid] = str(e)
                    self.html_writer.write('Error: %s</br>\n' % str(e))
                    logging.debug("C%05d: %s" % (cid, str(e)))
                    continue

                groupvector = decomposition.AsVector()
                try:
                    dG0 = self.groupvec2val(groupvector)
                except GroupMissingTrainDataError as e:
                    # in this case we do not care if a compound violated the group
                    # conservation laws because it might cancel out later when we 
                    # use it to calculate reactions.
                    dG0 = e.value
                    self.html_writer.write('Warning: %s</br>\n' % str(e))
                    logging.debug("C%05d: %s" % (cid, str(e)))
                    self.cid2error[cid] = str(e)

                source_string = "Group Contribution"
                
                if self.transformed:
                    diss_table.SetTransformedFormationEnergy(dG0_tag=dG0, 
                        pH=default_pH, I=default_I, pMg=default_pMg, T=default_T)
                else:
                    if nH != groupvector.Hydrogens() or nMg != groupvector.Magnesiums():
                        err_msg = "The most abundant pseudoisomer is [nH=%d, nMg=%d], " \
                            "but the decomposition has [nH=%d, nMg=%d]. Skipping..." \
                            "" % (nH, nMg, groupvector.Hydrogens(), groupvector.Magnesiums())
                        self.html_writer.write('ERROR: %s</br>\n' % err_msg)
                        self.cid2error[cid] = err_msg
                        continue
                    diss_table.SetFormationEnergyByNumHydrogens(dG0=dG0, nH=nH, nMg=nMg)
                pmap = diss_table.GetPseudoisomerMap()

            #if diss_table is not None:
            #    self.html_writer.write('<b>Dissociation table:</b></br>\n')
            #    diss_table.WriteToHTML(self.html_writer)
            
            if groupvector:
                self.html_writer.write('Group vector = %s</br>\n' % str(groupvector))
                self.html_writer.write('nH = %d, charge = %d, nMg = %d</br>\n' % 
                                       (groupvector.Hydrogens(),
                                        groupvector.NetCharge(),
                                        groupvector.Magnesiums()))

            self.SetPseudoisomerMap(cid, pmap)
            self.cid2groupvec[cid] = groupvector
            self.cid2source_string[cid] = source_string
            #self.html_writer.write('<b>Pseudoisomer table:</b></br>\n')
            #pmap.WriteToHTML(self.html_writer)
            
        self.html_writer.div_end()
        
        logging.info("Writing the results to the database")
        G.ToDatabase(self.db, table_name=self.PMAP_TABLE_NAME, 
                     error_table_name=self.ERROR_TABLE_NAME)
        self.db.CreateTable(self.GROUPVEC_TABLE_NAME, 'cid TEXT, groupvec TEXT')
        for cid, groupvec in self.cid2groupvec.iteritems():
            self.db.Insert(self.GROUPVEC_TABLE_NAME, [cid, groupvec.ToJSONString()])
        self.db.Commit()
        self.KeggErrorReport()

    def FromDatabase(self, db, table_name):
        # Read the formation energies of all pseudoisomers from the database
        reader = db.DictReader(table_name)
        PsuedoisomerTableThermodynamics._FromDictReader(
            reader, self, label=None, name="Group Contribution",
            warn_for_conflicting_refs=False)
        
        # Read the error messages from the database
        self.cid2error = {}            
        for row in db.DictReader(self.ERROR_TABLE_NAME):
            cid = int(row['cid'])
            if not cid:
                continue
            self.cid2error[cid] = row['error']

        # Read the group-vectors from the database
        self.cid2groupvec = {}
        for row in self.db.DictReader(self.GROUPVEC_TABLE_NAME):
            cid = int(row['cid'])
            groupvec = GroupVector.FromJSONString(self.groups_data, row['groupvec'])
            self.cid2groupvec[cid] = groupvec
    
    def Reaction2GroupVector(self, sparse):
        total_groupvec = GroupVector(self.groups_data)
        for cid, coeff in sparse.iteritems():
            if cid not in self.cid2groupvec: # one of the compounds has no observed or estimated formation energy
                raise MissingReactionEnergy("%s (C%05d) has no GroupVector" % 
                    (self.kegg.cid2name(cid), cid), sparse)
            groupvec = self.cid2groupvec[cid]
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
        nullspace_mat = []
        for row in self.db.DictReader(self.CONTRIBUTION_TABLE_NAME):
            self.group_contributions.append(row['dG0_gr'])
            nullspace_mat.append([float(x) for x in row['nullspace'].split(',')])
        self.group_contributions = pylab.array([self.group_contributions]).T
        self.group_nullspace = pylab.array(nullspace_mat).T
                        
    def GetGroupContribution(self, name, nH, z, nMg=0):
        gr = Group(None, name, nH, z, nMg)
        gid = self.groups_data.Index(gr)
        dG0_gr = self.group_contributions[gid, 0]
        return dG0_gr
    
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

    def KeggErrorReport(self):
        error_strings = ['kernel', 'decompose', 'explicit', 'provided', 'dissociation']
        query = ' UNION '.join(["SELECT '%s', COUNT(*) FROM %s WHERE error LIKE '%%%s%%'" % 
                                (e, self.ERROR_TABLE_NAME, e) for e in error_strings])
        self.db.Query2HTML(self.html_writer, query, ['Error', 'Count'])

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
            G.LoadGroupsFromFile()
            G.Train()
            G.WriteRegressionReport()
            G.analyze_training_set()
        else:
            G.init()

        if not options.train_only:
            G.EstimateKeggCids()

