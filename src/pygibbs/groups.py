#!/usr/bin/python

import csv
import logging
import types
import json

from copy import deepcopy
import pylab
import sys
from pygibbs.thermodynamic_constants import R, default_pH, default_T,\
    default_c0, dG0_f_Mg, correction_function
from pygibbs.thermodynamics import MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException,\
    KeggReactionNotBalancedException
from pygibbs.groups_data import Group, GroupsData
from pygibbs.pseudoisomer import PseudoisomerMap
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
from pygibbs.nist import Nist

class GroupMissingTrainDataError(Exception):
    
    def __init__(self, value, kernel_rows=[]):
        self.value = value
        self.kernel_rows = kernel_rows
    
    def __str__(self):
        if type(self.value) == types.StringType:
            return self.value
        else:
            return repr(self.value)
    
    def Explain(self, gc):
        missing_single_groups = []
        for i in self.kernel_rows:
            nonzero_columns = pylab.find(abs(gc.group_nullspace[i, :]) > 1e-10)
            if len(nonzero_columns) == 1:
                missing_single_groups.append(nonzero_columns[0])
            else:
                return 'contains group combinations that are not covered by ' + \
                    'training data'
                
        return 'contains missing groups: ' + ", ".join(
            [str(gc.groups_data.all_groups[j]) for j in missing_single_groups])
        

class GroupContribution(PsuedoisomerTableThermodynamics):    
    def __init__(self, db, html_writer=None, kegg=None):
        """Construct a GroupContribution instance.
        
        Args:
            db: the database handle to read from.
            html_writer: the HtmlWriter to write to.
            kegg: a Kegg instance if you don't want to use the default one.
        """
        PsuedoisomerTableThermodynamics.__init__(self, name="Group Contribution")
        self.db = db
        self.html_writer = html_writer or NullHtmlWriter()

        self.kegg = kegg or Kegg.getInstance()
        self.bounds = deepcopy(self.kegg.cid2bounds)
        self.sparse_kernel = False
        self.verify_kernel_orthogonality = True

        self.group_nullspace = None
        self.group_contributions = None
        self.obs_collection = None
        
        self.cid2error = None
        self.cid2groupvec = None
        
        self.PMAP_TABLE_NAME = 'gc_pseudoisomers'
        self.ERROR_TABLE_NAME = 'gc_errors'
        self.GROUPVEC_TABLE_NAME = 'gc_groupvector'
        self.NULLSPACE_TABLE_NAME = 'gc_nullspace'
        self.CONTRIBUTION_TABLE_NAME = 'gc_contribution'
                    
    def init(self):
        self.load_groups()
        self.LoadContributionsFromDB()
        self.load_concentrations()
        self.load_training_data()
        self.FromDatabase(self.db, table_name=self.PMAP_TABLE_NAME)
        
    def write_data_to_json(self, json_fname, kegg):
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
        
    def load_groups(self, group_fname=None):
        if group_fname:
            self.groups_data = GroupsData.FromGroupsFile(group_fname)
            self.groups_data.ToDatabase(self.db)
        else:
            self.groups_data = GroupsData.FromDatabase(self.db)
        self.group_decomposer = GroupDecomposer(self.groups_data)

    def read_training_data_reaction(self):
        self.html_writer.write('</br><b>List of NIST reactions for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.obs_collection.AddNistDatabase()
        self.html_writer.div_end()
        
    def read_training_data_formation(self):
        """
            Finds all the compounds which have a valid dG0 in the formation energies file,
            and adds the qualifying observations to obs_collection.
        """
        self.html_writer.write('<br><b>List of compounds for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.obs_collection.AddFormationEnergies()
        self.html_writer.div_end()
        
    def read_training_data_pKa(self):
        self.html_writer.write('<b>List of pKa for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.obs_collection.AddDissociationTable()
        self.html_writer.div_end()
            
    def load_training_data(self):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        self.obs_collection.FromDatabase()
        self.group_matrix, self.obs, self.obs_types, self.obs_names = self.obs_collection.GetRegressionData()

    def train(self, FromFiles=True):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        if FromFiles:
            #self.read_training_data_pKa()
            self.read_training_data_formation()
            self.read_training_data_reaction()
            self.obs_collection.ToDatabase()
            self.obs_collection.ToCSV('../res/observations.csv')
        else:
            self.obs_collection.FromDatabase()
            
        self.group_matrix, self.obs, self.obs_types, self.obs_names = self.obs_collection.GetRegressionData()
        
        self.group_contributions, self.group_nullspace = self.RunLinearRegression()
        self.SaveContributionsToDB()
            
    def RunLinearRegression(self):
        group_contributions, _nullspace = LinearRegression.LeastSquares(
                        self.group_matrix, self.obs, reduced_row_echlon=False)
        
        nullspace = _nullspace
        if self.sparse_kernel:
            try:
                nullspace = SparseKernel(self.group_matrix).Solve()
            except CplexNotInstalledError:
                logging.warning("CPLEX is not installed on this system, using a non-sparse"
                                " method for describing the Kernel of the group matrix")
            
        return list(group_contributions.flat), nullspace
    
    def does_table_exist(self, table_name):
        for unused_ in self.db.Execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def load_concentrations(self):
        self.media_list = []
        if (self.does_table_exist('compound_abundance')):
            for row in self.db.Execute("SELECT media FROM compound_abundance GROUP BY media"):
                self.media_list.append(row[0])
                
            self.cid2conc = {}
            for row in self.db.Execute("SELECT cid, media, concentration FROM compound_abundance"):
                cid, media, conc = row
                self.cid2conc[(cid, media)] = conc # in [M]

    def get_concentration(self, cid, c0=default_c0, media=None):
        if cid == 1: # the concentration of water must always be 1
            return 1
        if not media:
            return c0 # Standard conditions = 1 [M]
        return self.cid2conc.get((cid, media), c0)

    def get_concentration_list(self, cid, c0=default_c0):
        """
            return a list of pairs of media names and concentrations of the provided CID
        """
        c_list = []
        for media in self.media_list:
            if ((cid, media) in self.cid2conc):
                c_list.append((media, self.cid2conc[(cid, media)]))
        return c_list        
            
    def get_pseudoisomers(self, mol):
        pseudoisomers = set()
        decomposition = self.group_decomposer.Decompose(mol, ignore_protonations=True)
        
        for groupvec in decomposition.PseudoisomerVectors():
            nH = groupvec.Hydrogens()
            z = groupvec.NetCharge()
            mgs = groupvec.Magnesiums()
            pseudoisomers.add((nH, z, mgs))
        return sorted(list(pseudoisomers))
        
    def groupvec2val(self, groupvec):
        if self.group_contributions == None or self.group_nullspace == None:
            raise Exception("You need to first Train the system before using it to estimate values")

        if self.verify_kernel_orthogonality:
            v = abs(pylab.dot(self.group_nullspace, pylab.array(groupvec)))
            k_list = [i for i in pylab.find(v > 1e-10)]
            if k_list:
                raise GroupMissingTrainDataError("can't estimate because the input "
                    "is not orthogonal to the kernel", k_list)
        return pylab.dot(groupvec, self.group_contributions)
    
    def write_regression_report(self, include_raw_data=False, T=default_T, pH=default_pH):        
        if include_raw_data:
            for table_name in ['train_groups', 'group_observations']:
                self.html_writer.write('<b>%s</b>' % table_name)
                div_id = self.html_writer.insert_toggle()
                self.html_writer.div_start(div_id)
                self.html_writer.write('</br>')
                self.db.Table2HTML(self.html_writer, 'table_name')
                self.div_end()
            
        self.html_writer.write('</br><b>Regression report</b>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)
        self.html_writer.write_ul(['observations: %d' % self.group_matrix.shape[0],
           'groups: %d' % self.group_matrix.shape[1],
           'rank: %d' % LinearRegression.Rank(self.group_matrix)])
        self.html_writer.write('</br><table border="1">\n<tr>'
            '<td width="5%%">#</td>'
            '<td width="20%%">Name</td>'
            '<td width="5%%">&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td>'
            '<td width="70%%">Group Vector</td></tr>')
        for i in xrange(self.group_matrix.shape[0]):
            group_vector = self.group_matrix[i, :]
            nonzero_columns = pylab.find(abs(group_vector) > 1e-10)
            s_vector = " | ".join(["%g : %s" % (group_vector[0, j], self.groups_data.all_groups[j].name)
                for j in nonzero_columns])
            #s_vector = ['%s : %d' % (self.groups_data.all_groups[j].name, group_vector[0, j])
            #     for j in xrange(self.group_matrix.shape[1]) if group_vector[0, j] != 0]
            self.html_writer.write('<tr>' +
                '<td width="5%%">%d</td>' % i +
                '<td width="20%%">%s</td>' % self.obs_names[i] + 
                '<td width="5%%">%8.2f</td>' % self.obs[i] +
                '<td width="70%%">%s</td></tr>' % s_vector)
        self.html_writer.write('</table>')
        self.html_writer.div_end()
        
        self.html_writer.write('</br><b>Group Contributions</b>\n')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)
        self.html_writer.write('</br><font size="1">\n')
        dict_list = []
        for j, group in enumerate(self.groups_data.all_groups):
            dG0_gr = self.group_contributions[j]
            obs_lists_dict = {'acid-base':[], 'formation':[], 'reaction':[]}
            for k in pylab.find(self.group_matrix[:, j]):
                obs_lists_dict[self.obs_types[k]].append(self.obs_names[k])
            d = {"#":"%d" % j, "Group Name":group.name, 
                 "nH":group.hydrogens, "charge":group.charge, "nMg":group.nMg,
                 "&#x394;<sub>gr</sub>G [kJ/mol]":"%8.2f" % dG0_gr,
                 "dissociations":' | '.join(obs_lists_dict['acid-base']),
                 "formations":' | '.join(obs_lists_dict['formation']),
                 "reactions":' | '.join(obs_lists_dict['reaction'])}
            dict_list.append(d)
        self.html_writer.write_table(dict_list, headers=["#", "Group Name",
            "nH", "charge", "nMg", "&#x394;<sub>gr</sub>G [kJ/mol]", 
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
            row = {'name':self.obs_names[i], 'obs':self.obs[i]}
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.obs_types[i] not in ['formation', 'reaction']:
                continue
            
            row['fit'] = pylab.dot(self.group_matrix[i, :], self.group_contributions)[0, 0]
            row['fit_resid'] = row['fit']-row['obs'] 
            deviations.append(row)
            logging.info('Fit Error = %.1f' % (row['fit_resid']))

            # leave out the row corresponding with observation 'i'
            logging.info('Cross validation, leaving-one-out: ' + self.obs_names[i])
            subset = range(n_obs)
            subset.pop(i)
            loo_group_contributions, loo_nullspace = LinearRegression.LeastSquares(
                self.group_matrix[subset, :], self.obs[subset])
            
            if loo_nullspace.shape[0] > self.group_nullspace.shape[0]:
                logging.warning('example %d is not linearly dependent in the other examples' % i)
                continue
            row['loo'] = pylab.dot(self.group_matrix[i, :], loo_group_contributions)[0, 0]
            row['loo_resid'] = row['loo'] - row['obs']
            loo_deviations.append(row)
            logging.info('LOO Error = %.1f' % row['loo_resid'])
        
        logging.info("writing the table of estimation errors for each compound")
        self.html_writer.write('</br><b>Cross-validation table</b>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.div_start(div_id)

        obs_vec = pylab.array([row['obs'] for row in deviations])
        fit_vec = pylab.array([row['fit'] for row in deviations])
        
        resid_vec = pylab.array([row['fit_resid'] for row in deviations])
        rmse = pylab.rms_flat(resid_vec)
        
        loo_resid_vec = pylab.array([row['loo_resid'] for row in loo_deviations])
        loo_rmse = pylab.rms_flat(loo_resid_vec)

        self.html_writer.write_ul(['fit_rmse(pred) = %.2f kJ/mol' % rmse,
                                   'loo_rmse(pred) = %.2f kJ/mol' % loo_rmse])
        logging.info("Goodness of fit: RMSE = %.2f kJ/mol" % rmse)
        logging.info("Leave-one-out test: RMSE = %.2f kJ/mol" % loo_rmse)

        self.html_writer.write('<font size="1">\n')
        self.html_writer.table_start()
        headers = ['Compound Name',
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

    def Mol2PseudoisomerMap(self, mol, ignore_protonations=False):
        decomposition = self.Mol2Decomposition(mol, ignore_protonations)
        return self.GroupDecomposition2PseudoisomerMap(decomposition)
    
    def GroupDecomposition2PseudoisomerMap(self, decomposition):
        all_groupvecs = decomposition.PseudoisomerVectors()
        if not all_groupvecs:
            raise GroupDecompositionError('Found no pseudoisomers for %s'
                                          % str(decomposition.mol))

        kernel_rows = set()
        pmap = PseudoisomerMap()
        for groupvec in all_groupvecs:
            try:
                dG0 = self.groupvec2val(groupvec)
                pmap.AddGroupVector(groupvec, dG0)
            except GroupMissingTrainDataError as e:
                kernel_rows.update(e.kernel_rows)
        if pmap.Empty():
            raise GroupMissingTrainDataError("All species of %s have missing groups: " % 
                                             decomposition.mol, kernel_rows)            

        pmap.Squeeze()
        return pmap

    def GroupVectorToTransformedGibbsEnergy(self, groupvector,
            pH=None, pMg=None, I=None, T=None):
        pH = pH or self.pH
        pMg = pMg or self.pMg
        I = I or self.I
        T = T or self.T

        nH = groupvector.Hydrogens()
        z = groupvector.NetCharge()
        nMg = groupvector.Magnesiums()
        dG0 = self.groupvec2val(groupvector)
        dG0_tag = dG0 + correction_function(nH, z, nMg, pH, pMg, I, T)
        return dG0_tag
        
    def Decomposition2Psuedoisomers(self, decomposition, ignore_protonations=True):
        """
            Given a decomposition, returns a list of possible pseudoisomers
            together with their relative abundance in the conditions specified.
            Also checks that all pseudoisomers have a legal nH and charge,
            i.e. nH - z = const.
            
            Returns:
                A pair (pmap, group_vector), where the pmap is for all the
                pseudoisomers possible, and the group_vector is for the most
                abundant pseudoisomer at the specified conditions.
        """
        nH = decomposition.Hydrogens()
        z = decomposition.NetCharge()
        nMg = decomposition.Magnesiums()
        
        all_groupvecs = decomposition.PseudoisomerVectors()
        if not all_groupvecs:
            raise GroupDecompositionError('Found no pseudoisomers for %s'
                                          % str(decomposition.mol))

        kernel_rows = set()
        nH_and_nMg_to_species = {}
        for groupvector in all_groupvecs:
            try:
                d = {}
                d['group vector'] = groupvector
                d['nH'] = groupvector.Hydrogens()
                d['charge'] = groupvector.NetCharge()
                d['nMg'] = groupvector.Magnesiums()
                if (d['nH'] + 2*d['nMg'] - d['charge']) != (nH + 2*nMg - z):
                    continue
                
                d['dG0'] = self.groupvec2val(groupvector)
                if (d['nH'], d['nMg']) not in nH_and_nMg_to_species or \
                        d['dG0'] < nH_and_nMg_to_species[d['nH'], d['nMg']]['dG0']:
                    nH_and_nMg_to_species[d['nH'], d['nMg']] = d
                
            except GroupMissingTrainDataError as e:
                kernel_rows.update(e.kernel_rows)

        if not nH_and_nMg_to_species:
            raise GroupMissingTrainDataError("All species of %s have missing groups: " % 
                                             decomposition.mol, kernel_rows)
        
        pmap = PseudoisomerMap()
        for d in nH_and_nMg_to_species.values():
            pmap.Add(d['nH'], d['charge'], d['nMg'], d['dG0'], 
                     ref='Group Contribution')
        pmap.FilterImprobablePseudoisomers(threshold=0.001, T=self.T)
        return pmap, nH_and_nMg_to_species
    
    def EstimateKeggCids(self, cid_list=None):
        """
            Uses the Group Contributions to estimate the entire set of compounds in KEGG,
            and then writes the results to the database as 'gc_pseudoisomers' table
            
            Options:
                override_all_observed_compounds - If True, any observed formation energy is 
                    used instead of the GC estimation. If False, only 'test' compounds are used.
        """
        logging.info("Estimating formation energies for all KEGG. Please be patient for a few minutes...")
        
        obs_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='testing')
            
        dissociation = self.dissociation = DissociationConstants.FromDatabase(
                                    self.db, 'dissociation_constants_chemaxon')
        
        self.cid2pmap_dict = {}
        self.cid2source_string = {}
        self.cid2error = {}
        self.cid2groupvec = {}
        
        for cid in sorted(cid_list or self.kegg.get_all_cids()):
            self.html_writer.write('<p>\n')
            self.html_writer.write('C%05d - %s\n' % (cid, self.kegg.cid2name(cid)))

            mol = None
            diss_table = None
            if cid in obs_species.get_all_cids():
                pmap = obs_species.cid2PseudoisomerMap(cid)
                pmatrix = pmap.ToMatrix() # returns a list of (nH, z, nMg, dG0)
                if len(pmatrix) != 1:
                    raise Exception("C%05d has multiple training species" % cid)
                nH, z, nMg, dG0 = pmatrix[0]
                groupvector = GroupVector(self.groups_data) # use an empty group vector
                source_string = obs_species.cid2SourceString(cid)
            else:
                diss_table = dissociation.GetDissociationTable(cid, 
                                                               create_if_missing=False)
                if diss_table is None:
                    continue # TODO: replace this with the method that uses ChamAxon to find pKas
                mol = diss_table.GetMostAbundantMol(pH=default_pH, I=0, 
                                                    pMg=14, T=default_T)
                if mol is None:
                    continue
                try:
                    decomposition = self.Mol2Decomposition(mol,
                                                           ignore_protonations=False)
                except (GroupDecompositionError, GroupMissingTrainDataError) as e:
                    self.cid2error[cid] = str(e)
                    self.html_writer.write('</br>%s\n' % str(e))
                    continue

                groupvector = decomposition.AsVector()
                nH = groupvector.Hydrogens()
                z = groupvector.NetCharge()
                nMg = groupvector.Magnesiums()
                dG0 = self.groupvec2val(groupvector)
                source_string = "Group Contribution"

            self.html_writer.insert_toggle('C%05d' % cid)
            self.html_writer.div_start('C%05d' % cid)
            if mol:
                self.html_writer.write(mol.ToSVG() + '</br>\n')
            if groupvector:
                self.html_writer.write('Group vector = %s</br>\n' % str(groupvector))
                self.html_writer.write('nH = %d, charge = %d, nMg = %d</br>\n' % 
                                       (groupvector.Hydrogens(),
                                        groupvector.NetCharge(),
                                        groupvector.Magnesiums()))
            if diss_table:
                self.html_writer.write('<h2>ChemAxon dissociation constants:</h2>\n')
                diss_table.WriteToHTML(self.html_writer)
                self.html_writer.write('</br>\n')
            
                diss_table.SetFormationEnergyByNumHydrogens(dG0=dG0, nH=nH, nMg=nMg)
                pmap = diss_table.GetPseudoisomerMap()
            else:
                pmap = PseudoisomerMap(nH=nH, z=z, nMg=nMg, dG0=dG0,
                                       ref=source_string)
            self.SetPseudoisomerMap(cid, pmap)
            self.cid2groupvec[cid] = groupvector
            self.cid2source_string[cid] = source_string

            self.html_writer.write('<h2>Pseudoisomer table:</h2>\n')
            pmap.WriteToHTML(self.html_writer)
            
            self.html_writer.div_end()
            self.html_writer.write('</p>\n')
        
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
                raise MissingReactionEnergy("C%05d has no GroupVector" % cid, sparse)
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
        elif cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy(self.cid2error.get(cid, ""), cid)
        
    def SaveContributionsToDB(self):
        logging.info("storing the group contribution data in the database")
        self.db.CreateTable(self.CONTRIBUTION_TABLE_NAME,
                            'gid INT, name TEXT, protons INT, charge INT, '
                            'nMg INT, dG0_gr REAL, nullspace TEXT')
        
        for j, dG0_gr in enumerate(self.group_contributions):
            group = self.groups_data.all_groups[j]
            nullspace_str = ','.join(["%.2f" % x for x in self.group_nullspace[:, j]])
            self.db.Insert(self.CONTRIBUTION_TABLE_NAME, [j, group.name, group.hydrogens,
                                            group.charge, group.nMg, dG0_gr,
                                            nullspace_str])
            
        self.db.CreateTable(self.NULLSPACE_TABLE_NAME, 'dimension INT, group_vector TEXT')
        for i in xrange(self.group_nullspace.shape[0]):
            groupvec = GroupVector(self.groups_data, self.group_nullspace[i, :])
            self.db.Insert(self.NULLSPACE_TABLE_NAME, 
                           [i, groupvec.ToJSONString()])

        self.db.Commit()

    def GroupMatrixRowToString(self, r):
        nonzero_columns = pylab.find(abs(r) > 1e-10)
        return " | ".join(["%g : %s" % (r[j], self.groups_data.all_groups[j].name)
                           for j in nonzero_columns])
            
    def LoadContributionsFromDB(self):
        logging.info("loading the group contribution data from the database")
        self.group_contributions = []
        nullspace_mat = []
        for row in self.db.Execute("SELECT * FROM gc_contribution ORDER BY gid"):
            (unused_gid, unused_name, unused_protons, unused_charge, unused_mg, dG0_gr, nullspace) = row
            self.group_contributions.append(dG0_gr)
            nullspace_mat.append([float(x) for x in nullspace.split(',')])
        self.group_nullspace = pylab.matrix(nullspace_mat).T
                        
    def read_compound_abundance(self, filename):
        self.db.CreateTable('compound_abundance', 'cid INT, media TEXT, concentration REAL')
        for row in csv.DictReader(open(filename, 'r')):
            if not row['Kegg ID']:
                continue
            cid = int(row['Kegg ID'])
            if row['use'] != '1':
                continue
            try:
                self.db.Insert('compound_abundance', [cid, "glucose", float(row['Glucose'])])
            except ValueError:
                pass
            try:
                self.db.Insert('compound_abundance', [cid, "glycerol", float(row['Glycerol'])])
            except ValueError:
                pass
            try:
                self.db.Insert('compound_abundance', [cid, "acetate", float(row['Acetate'])])
            except ValueError:
                pass
        self.db.Commit()
        self.load_concentrations()

    def analyze_decomposition_cid(self, cid):
        return self.analyze_decomposition(self.kegg.cid2mol(cid))

    def analyze_decomposition(self, mol, ignore_protonations=False, strict=False):
        return self.group_decomposer.Decompose(mol, ignore_protonations, strict).ToTableString()
    
    def get_group_contribution(self, name, nH, z, nMg=0):
        gr = Group(None, name, nH, z, nMg)
        gid = self.groups_data.Index(gr)
        dG0_gr = self.group_contributions[gid]
        return dG0_gr
    
    def GetpKa_group(self, name, nH, z, nMg=0):
        dG0_gr_deprotonated = self.get_group_contribution(name, nH-1, z-1, nMg)
        dG0_gr_protonated = self.get_group_contribution(name, nH, z, nMg)
        if not dG0_gr_deprotonated or not dG0_gr_protonated:
            return None
        pKa = (dG0_gr_deprotonated - dG0_gr_protonated)/(R*default_T*pylab.log(10))
        return pKa
    
    def GetpK_Mg_group(self, name, nH, z, nMg):
        dG0_gr_without_Mg = self.get_group_contribution(name, nH, z-2, nMg-1)
        dG0_gr_with_Mg = self.get_group_contribution(name, nH, z, nMg)
        if not dG0_gr_without_Mg or not dG0_gr_with_Mg:
            return None
        pK_Mg = (dG0_gr_without_Mg + dG0_f_Mg - dG0_gr_with_Mg)/(R*default_T*pylab.log(10))
        return pK_Mg

    def KeggErrorReport(self):
        error_strings = ['kernel', 'decompose', 'explicit', 'provided', 'dissociation']
        query = ' union '.join(["select '" + e + 
            "', count(*) from gc_errors where error like '%%" + 
            e + "%%'" for e in error_strings])
        self.db.Query2HTML(self.html_writer, query, ['Error', 'Count'])

    def AnalyzeSingleCompound(self, mol, ignore_protonations=False):
        print 'Analyzing %s:' % str(mol)
        try:
            decomposition = self.Mol2Decomposition(mol, ignore_protonations=ignore_protonations)
            print decomposition.ToTableString()
            print 'nH =', decomposition.Hydrogens()
            print 'z =', decomposition.NetCharge()
            print 'nMg = ', decomposition.Magnesiums()
        except GroupDecompositionError as e:
            print "Cannot decompose compound to groups: " + str(e)
        mol.Draw()

    def EstimateCompoundDissociation(self, mol):
        decomposition = self.Mol2Decomposition(mol, ignore_protonations=True)
        groupvector = self.GroupDecomposition2MostAbuntantPseudoisomer(decomposition, pH=7)
        pmap = self.GroupDecomposition2PseudoisomerMap(decomposition)

        nH = groupvector.Hydrogens()
        z = groupvector.NetCharge()
        nMg = groupvector.Magnesiums()
        dG0 = self.groupvec2val(groupvector)
        
        print "- Most Abundante: nH=%d, z=%d, nMg=%d, dG0=%.1f" % (nH, z, nMg, dG0)
        print "- PsuedoisomerMap:\n", pmap
        print "- Dissociation Table:"
        for nH_below, nH_above, pKa in pmap.GetAllpKas():
            print "[nH=%d -> nH=%d] : pKa = %.1f" % (nH_below, nH_above, pKa)

    def EstimateGroupDissociation(self):
        name2group_list = {}
        for i, group in enumerate(self.groups_data.groups):
            dG0 = self.group_contributions[i]
            name2group_list.setdefault(group.name, []).append((group, dG0))
        
        for _name, group_list in name2group_list.iteritems():
            if len(group_list) > 1:
                print '-' * 50
                for group, dG0 in group_list:
                    print group, dG0

#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
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
    opt_parser.add_option("-p", "--pka", action="store_true",
                          dest="pka", default=False,
                          help="A flag for running a pKa analysis of the groups")
    return opt_parser

if __name__ == '__main__':
    options, _ = MakeOpts().parse_args(sys.argv)
    util._mkdir('../res')
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    if options.pka: # use flag -p of --pka
        G = GroupContribution(db=db)
        G.init()
        G.EstimateGroupDissociation()
    elif options.smiles: # -s <SMILES>
        G = GroupContribution(db=db)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        print 'Analyzing SMILES %s:' % (options.smiles)
        mol = Molecule.FromSmiles(options.smiles)
        G.AnalyzeSingleCompound(mol, ignore_protonations=False)
    elif options.inchi: # -i <INCHI>
        G = GroupContribution(db=db)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        print 'Analyzing InChI %s:' % (options.inchi)
        mol = Molecule.FromInChI(options.inchi)
        G.AnalyzeSingleCompound(mol, ignore_protonations=False)
    elif options.cid: # -c <CID>
        G = GroupContribution(db=db)
        G.init()
        mol = G.kegg.cid2mol(options.cid)
        G.EstimateCompoundDissociation(mol)
        G.AnalyzeSingleCompound(mol, ignore_protonations=True)
    elif options.rid: # -r <RID>
        G = GroupContribution(db=db)
        G.init()
        reaction = G.kegg.rid2reaction(options.rid)
        dG0_r = reaction.PredictReactionEnergy(G)
        print "R%05d = %.2f" % (options.rid, dG0_r)
    else:
        # use the flag -i or --train for train only
        # use the flag -e or --test for test only
        
        if options.test_only:
            html_writer = HtmlWriter('../res/groups_test.html')
        elif options.train_only:
            html_writer = HtmlWriter('../res/groups_train.html')
        else:
            html_writer = HtmlWriter('../res/groups.html')
            
        G = GroupContribution(db=db, html_writer=html_writer)
        G.verify_kernel_orthogonality = False
        if not options.test_only:
            G.load_groups("../data/thermodynamics/groups_species.csv")
            G.train()
            G.write_regression_report()
            G.analyze_training_set()
        else:
            G.init()

        if not options.train_only:
            nist = Nist()
            G.EstimateKeggCids(nist.GetAllCids())
            #G.verify_kernel_orthogonality = True
            #G.EstimateKeggRids()

