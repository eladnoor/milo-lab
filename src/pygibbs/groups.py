#!/usr/bin/python

import csv
import logging
import types
import json

from copy import deepcopy
import pylab
import sys
from pygibbs.thermodynamic_constants import R, default_pH, default_T, default_c0, dG0_f_Mg,\
    default_I
from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException
from pygibbs.groups_data import Group, GroupsData
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs import templates
from toolbox import util
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.sparse_kernel import SparseKernel, CplexNotInstalledError
from toolbox.molecule import Molecule
from pygibbs.group_observation import GroupObervationCollection
from pygibbs.dissociation_constants import DissociationConstants,\
    MissingDissociationConstantError

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
        

class GroupContribution(Thermodynamics):    
    def __init__(self, db, html_writer=None, kegg=None):
        """Construct a GroupContribution instance.
        
        Args:
            db: the database handle to read from.
            html_writer: the HtmlWriter to write to.
            kegg: a Kegg instance if you don't want to use the default one.
        """
        Thermodynamics.__init__(self)
        self.db = db
        self.html_writer = html_writer or NullHtmlWriter()

        self.kegg = kegg or Kegg.getInstance()
        self.bounds = deepcopy(self.kegg.cid2bounds)
        self.sparse_kernel = True
        self.verify_kernel_orthogonality = True

        self.group_nullspace = None
        self.group_contributions = None
        self.obs_collection = None
        
        self.cid2pmap = None
        self.cid2error = None
            
    def write_gc_tables(self):
        table_names = ["groups", "contribution", "observation"]
        self.html_writer.write('<ul>\n')
        for table_name in table_names:
            self.html_writer.write('  <li><a href="#%s">Table %s from the database</a></li>\n' % (table_name, table_name))
        self.html_writer.write('</ul>\n')
        for table_name in table_names:
            self.html_writer.write('<p><h2><a name="%s">Table %s:</a></h2>\n' % (table_name, table_name))
            self.db.Table2HTML(self.html_writer, table_name)
            self.html_writer.write('</p>\n')
    
    def init(self):
        self.load_groups()
        self.LoadContributionsFromDB()
        self.load_concentrations()
        self.load_training_data()
        
    def save_cid2pmap(self):
        logging.info('Calculating the table of chemical formation energies for all KEGG compounds')
        self.db.CreateTable('gc_cid2prm', 'cid INT, nH INT, z INT, nMg INT, dG0 REAL, estimated BOOL')
        self.db.CreateTable('gc_cid2error', 'cid INT, error TEXT')

        compounds = []

        for cid in self.kegg.get_all_cids():
            cdict = {'cid': cid, 'measured_pmap': None,
                     'estimated_pmap': None, 'compound': None}
            
            if cid % 100 == 1:
                logging.info('Saving KEGG Compound C%05d' % cid)
            
            # If the compound is measured:
            if cid in self.obs_collection.cid2pmap:
                pmap = self.obs_collection.cid2pmap[cid]
                cdict['measured_pmap'] = pmap
                for nH, z, nMg, dG0 in pmap.ToMatrix():
                    self.db.Insert('gc_cid2prm', [cid, nH, z, nMg, dG0, False])

            # Try to also estimate the dG0_f using Group Contribution:
            comp = self.kegg.cid2compound(cid)
            cdict['compound'] = comp
            error_str = None
            if not comp.inchi:
                error_str = 'no InChI exists'
            else:
                try:
                    decomposition = self.Mol2Decomposition(comp.GetMolecule(),
                                                           ignore_protonations=True)
                    cdict['decomposition'] = decomposition
                    pmap = self.GroupDecomposition2PseudoisomerMap(decomposition)
                    cdict['estimated_pmap'] = pmap
                    for nH, z, nMg, dG0 in pmap.ToMatrix():
                        self.db.Insert('gc_cid2prm', [cid, int(nH), int(z), int(nMg), dG0, True])
                except KeggParseException:
                    error_str = 'cannot determine molecular structure'
                except GroupDecompositionError:
                    error_str = 'cannot decompose into groups'
                except GroupMissingTrainDataError as e:
                    error_str = e.Explain(self)
            compounds.append(cdict)
            self.db.Insert('gc_cid2error', [cid, error_str])
        
        logging.info('Rendering the kegg_pmaps.html file')
        templates.render_to_file('kegg_pmaps.html', {'compounds': compounds},
                                 '../res/kegg_pmaps.html')
        logging.info('Writing the formation energies to the database')
        self.db.Commit()
        logging.info('DONE!')

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
        self.html_writer.write('<h2><a name=compounds>List of NIST reactions for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.obs_collection.AddNistDatabase()
        self.html_writer.write('</div>')
        
    def read_training_data_formation(self):
        """
            Finds all the compounds which have a valid dG0 in the formation energies file,
            and adds the qualifying observations to obs_collection.
        """
        self.html_writer.write('<h2><a name=compounds>List of compounds for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.obs_collection.AddFormationEnergies()
        self.html_writer.write('</div>')
        
    def read_training_data_pKa(self):
        self.html_writer.write('<h2><a name=compounds>List of pKa for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.obs_collection.AddDissociationTable()
        self.html_writer.write('</div>')
            
    def load_training_data(self):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        self.obs_collection.FromDatabase()
        self.group_matrix, self.obs, self.obs_types, self.obs_names = self.obs_collection.GetRegressionData()

    def train(self, FromFiles=True):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        if FromFiles:
            self.read_training_data_formation()
            self.read_training_data_pKa()
            self.read_training_data_reaction()
            self.obs_collection.ToDatabase()
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
        self.html_writer.write('<h2><a name=compounds>Regression report</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)

        self.html_writer.write_ul(['%d compounds' % self.group_matrix.shape[0],
                            '%d groups' % self.group_matrix.shape[1]])
        
        if include_raw_data:
            self.db.Table2HTML(self.html_writer, 'train_groups')
            self.db.Table2HTML(self.html_writer, 'group_observations')
            
        self.html_writer.write('<table border="1">\n<tr>'
                        '<td>Name</td>'
                        '<td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td>'
                        '<td>Group Vector</td></tr>')
        for i in xrange(self.group_matrix.shape[0]):
            self.html_writer.write('<tr><td>%s</td><td>%.2f</td><td>%s</td></tr>' % \
                            (self.obs_names[i], self.obs[i], 
                             ' '.join(['%d' % self.group_matrix[i, j] for j in xrange(self.group_matrix.shape[1])])))
        self.html_writer.write('</table>')
        
        self.html_writer.write('<h2><a name="group_contrib">Group Contributions</a></h2>\n')
        self.html_writer.write('<table border="1">')
        self.html_writer.write('  <tr><td>#</td><td>Group Name</td><td>nH</td><td>charge</td><td>nMg</td><td>&#x394;<sub>gr</sub>G [kJ/mol]</td><td>&#x394;<sub>gr</sub>G\' [kJ/mol]</td><td>Appears in compounds</td></tr>\n')
        for i, group in enumerate(self.groups_data.all_groups):
            dG0_gr = self.group_contributions[i]
            dG0_gr_tag = dG0_gr + R*T*pylab.log(10)*pH*group.hydrogens
            compound_list_str = ' | '.join([self.obs_names[k] for k in pylab.find(self.group_matrix[:, i])])
            self.html_writer.write('  <tr><td>%d</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%8.2f</td><td>%8.2f</td><td>%s</td></tr>\n' %
                            (i, group.name, group.hydrogens,
                             group.charge, group.nMg,
                             dG0_gr, dG0_gr_tag, compound_list_str))
        self.html_writer.write('</table>\n')
        self.html_writer.write('</div>\n')

        self.html_writer.write('<h2><a name="nullspace">Nullspace of regression matrix</a></h2>\n')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('<div id="%s" style="display:none">\n' % div_id)
        self.db.Table2HTML(self.html_writer, 'gc_nullspace')
        self.html_writer.write('</div>\n')

    def analyze_training_set(self):
        n_obs = self.group_matrix.shape[0]
        val_names = []
        val_obs = []
        val_est = []
        val_err = []
        deviations = []
        
        for i in xrange(n_obs):
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.obs_types[i] not in ['formation', 'reaction']:
                continue
            logging.info('Cross validation, leaving-one-out: ' + self.obs_names[i])
            
            # leave out the row corresponding with observation 'i'
            subset = range(n_obs)
            subset.pop(i)
            
            group_contributions, nullspace = LinearRegression.LeastSquares(
                self.group_matrix[subset, :], self.obs[subset])
            
            if nullspace.shape[0] > self.group_nullspace.shape[0]:
                logging.warning('# example %d is not linearly dependent in the other examples' % i)
                deviations.append((None, self.obs_names[i], "%.1f" % self.obs[i], "-", "-", 
                                   "linearly independent example", "-"))
                continue
            
            orig_estimation = pylab.dot(self.group_matrix[i, :], self.group_contributions)
            estimation = pylab.dot(self.group_matrix[i, :], group_contributions)[0, 0]
            error = self.obs[i] - estimation
            
            val_names.append(self.obs_names[i])
            val_obs.append(self.obs[i])
            val_est.append(estimation)
            val_err.append(error)
            logging.info('Error = %.1f' % error)
            deviations.append((abs(error), self.obs_names[i], 
                               "%.1f" % self.obs[i], 
                               "%.1f" % estimation,
                               "%.1f" % error, "",
                               "%.1f" % orig_estimation))
        
        logging.info("writing the table of estimation errors for each compound")
        self.html_writer.write('<h2><a name=compounds>Cross Validation Table</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)

        rmse = util.calc_rmse(val_obs, val_est)
        self.html_writer.write('<div><b>rmse(pred) = %.2f kJ/mol</b></div>\n' % rmse)
        logging.info("The RMSE of the predictions is %.2f kJ/mol" % rmse)

        self.html_writer.write('<table border="1">')
        self.html_writer.write('  <tr><td>Compound Name</td>'
                        '<td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td>'
                        '<td>&#x394;<sub>f</sub>G<sub>est</sub> [kJ/mol]</td>'
                        '<td>Error [kJ/mol]</td>'
                        '<td>Remark</td><td>&#x394;<sub>f</sub>G<sub>orig est</sub> [kJ/mol]</td></tr>\n')
        deviations.sort(reverse=True)
        for _, name, obs, est, err, remark, orig in deviations:
            self.html_writer.write('  <tr><td>' + '</td><td>'.join([name, obs, est, err, remark, orig]) + '</td></tr>\n')
        self.html_writer.write('</table>\n')
        self.html_writer.write('</div>\n')
        
        logging.info("Plotting graphs for observed vs. estimated")
        self.html_writer.write('<h2><a name=compounds>Cross Validation Figure 1</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)

        obs_vs_est_fig = pylab.figure()
        pylab.plot(val_obs, val_est, '.')
        pylab.xlabel('Observed (obs)')
        pylab.ylabel('Estimated (est)')
        pylab.hold(True)
        for i in xrange(len(val_names)):
            pylab.text(val_obs[i], val_est[i], val_names[i], fontsize=4)
        self.html_writer.write('<h3><a name="obs_vs_est">Observed vs. Estimated</a></h3>\n')
        self.html_writer.embed_matplotlib_figure(obs_vs_est_fig, width=1000, height=800)

        self.html_writer.write('</div>\n')
        self.html_writer.write('<h2><a name=compounds>Cross Validation Figure 2</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        
        obs_vs_err_fig = pylab.figure()
        pylab.plot(val_obs, val_err, '+')
        pylab.xlabel('Observed (obs)')
        pylab.ylabel('Estimation error (est - obs)')
        pylab.hold(True)
        for i in xrange(len(val_names)):
            pylab.text(val_obs[i], val_err[i], val_names[i], fontsize=4)
        self.html_writer.write('<h3><a name="obs_vs_err">Observed vs. Error</a></h3>\n')
        self.html_writer.embed_matplotlib_figure(obs_vs_err_fig, width=1000, height=800)
        self.html_writer.write('</div>\n')

    def get_all_cids(self):
        return sorted(self.kegg.get_all_cids())

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

    def EstimateKeggCids(self, override_all_observed_compounds=False):
        """
            Uses the Group Contributions to estimate the entire set of compounds in KEGG,
            and then writes the results to the database as 'gc_pseudoisomers' table
            
            Options:
                override_all_observed_compounds - If True, any observed formation energy is 
                    used instead of the GC estimation. If False, only 'test' compounds are used.
        """
        logging.info("Estimating formation energies for all KEGG. Please be patient for a few minutes...")
        
        # override the data of the 'test' set
        if override_all_observed_compounds:
            obs_species = PsuedoisomerTableThermodynamics.FromCsvFile(
                '../data/thermodynamics/formation_energies.csv')
        else:
            obs_species = PsuedoisomerTableThermodynamics.FromCsvFile(
                '../data/thermodynamics/formation_energies.csv', label='test')
            
        obs_cids = obs_species.get_all_cids()

        dissociation = DissociationConstants.FromFile()
        self.cid2pmap = {}
        self.cid2error = {}
        for cid in sorted(self.kegg.get_all_cids()):
            try:
                diss = dissociation.GetDissociationTable(cid, create_if_missing=False)
                if diss is None:
                    mol = self.kegg.cid2mol(cid)
                    decomposition = self.Mol2Decomposition(mol, ignore_protonations=True)
                    nH = decomposition.Hydrogens()
                    z = decomposition.NetCharge()
                    nMg = decomposition.Magnesiums()
                    dG0 = self.groupvec2val(decomposition.AsVector())
                    self.cid2pmap[cid] = PseudoisomerMap(nH=nH, z=z, nMg=nMg, dG0=dG0)
                    self.cid2source_string[cid] = "Group Contribution"
                else:
                    if cid in obs_cids:
                        pmap = obs_species.cid2PseudoisomerMap(cid)
                        pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
                        if len(pmatrix) != 1:
                            raise Exception("C%05d has more than one species in the training set" % cid)
                        nH, _z, nMg, dG0 = pmatrix[0]
                        self.cid2source_string[cid] = obs_species.cid2SourceString(cid)
                    else:
                        nH, nMg = diss.GetMostAbundantPseudoisomer(
                            pH=default_pH, I=default_I, pMg=14, T=default_T)
                        smiles = diss.GetSmiles(nH=nH, nMg=nMg)
                        if smiles:
                            mol = Molecule.FromSmiles(smiles)
                            decomposition = self.Mol2Decomposition(mol, ignore_protonations=False)
                            if nH != decomposition.Hydrogens():
                                raise GroupDecompositionError("The nH of the decomposition (%d) does not match the provided nH (%d)"
                                    % (decomposition.Hydrogens(), nH))
                        else:
                            mol = self.kegg.cid2mol(cid)
                            decomposition = self.Mol2Decomposition(mol, ignore_protonations=True)
                            nH = decomposition.Hydrogens()
                            nMg = decomposition.Magnesiums()
                        dG0 = self.groupvec2val(decomposition.AsVector())
                        self.cid2source_string[cid] = "Group Contribution"
                
                    diss.SetFormationEnergyByNumHydrogens(dG0=dG0, nH=nH, nMg=nMg)
                    self.cid2pmap[cid] = diss.GetPseudoisomerMap()
    
            except (KeggParseException, GroupDecompositionError, 
                    GroupMissingTrainDataError, MissingDissociationConstantError) as e:
                self.cid2error[cid] = str(e)
                continue
        
        logging.info("Writing the results to the database")
        G.ToDatabase(db, table_name='gc_pseudoisomers', error_table_name='gc_errors')
        self.KeggErrorReport()
            
    def cid2PseudoisomerMap(self, cid):
        """
            Returns a PseudoisomerMap according to a given CID.
        """
        if self.cid2pmap is None:
            raise Exception("You must run the method EstimateKeggCids before using this one.")
        if cid in self.cid2pmap:
            return self.cid2pmap[cid]
        else:
            raise MissingCompoundFormationEnergy(self.cid2error[cid])
        
    def SaveContributionsToDB(self):
        logging.info("storing the group contribution data in the database")
        self.db.CreateTable('gc_contribution', 'gid INT, name TEXT, protons INT, charge INT, nMg INT, dG0_gr REAL, nullspace TEXT')
        
        for j, dG0_gr in enumerate(self.group_contributions):
            group = self.groups_data.all_groups[j]
            nullspace_str = ','.join(["%.2f" % x for x in self.group_nullspace[:, j]])
            self.db.Insert('gc_contribution', [j, group.name, group.hydrogens,
                                            group.charge, group.nMg, dG0_gr,
                                            nullspace_str])
            
        self.db.CreateTable('gc_nullspace', 'dimension INT, group_vector TEXT')
        for i in xrange(self.group_nullspace.shape[0]):
                self.db.Insert('gc_nullspace', [i, self.NullspaceRowToString(i)])

        self.db.Commit()

    def NullspaceRowToString(self, i):
        nonzero_columns = pylab.find(abs(self.group_nullspace[i, :]) > 1e-10)
        return ",".join(["%g x %s" % (self.group_nullspace[i, j], 
            str(self.groups_data.all_groups[j])) for j in nonzero_columns])
            
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
        
#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
    
if __name__ == '__main__':
    db = SqliteDatabase('../res/gibbs.sqlite')
    if len(sys.argv) < 2:
        html_writer = HtmlWriter('../res/groups.html')
        G = GroupContribution(db=db, html_writer=html_writer)
        G.verify_kernel_orthogonality = False
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.train()
        G.write_regression_report()
        G.analyze_training_set()
        G.EstimateKeggCids()
    elif sys.argv[1] == 'test':
        html_writer = HtmlWriter('../res/groups.html')
        G = GroupContribution(db=db, html_writer=html_writer)
        G.init()
        G.EstimateKeggCids()
    else:
        G = GroupContribution(db=db)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        #G.init()

        mols = {}
        try:
            cid = int(sys.argv[1])
            mols[G.kegg.cid2smiles(cid)] = G.kegg.cid2mol(cid)
        except ValueError:
            mols['smiles'] = Molecule.FromSmiles(sys.argv[1])
        
        #mols['ATP'] = Molecule.FromSmiles('C(C1C(C(C(n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O')
        #mols['Tryptophan'] = Molecule.FromSmiles("c1ccc2c(c1)c(CC(C(=O)O)[NH3+])c[nH]2")
        #mols['Adenine'] = Molecule.FromSmiles('c1nc2c([NH2])[n]c[n-]c2n1')
        #mols['Glutamate'] = Molecule.FromSmiles('C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OP(=O)([O-])[O-])O)O)O)O')
        #mols['Acetyl-CoA [nH=36]'] = Molecule.FromSmiles('CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])(O)=O)n2cnc3c(N)[nH+]cnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')
        #mols['Acetyl-CoA [nH=35]'] = Molecule.FromSmiles('CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])(O)=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')
        #mols['Acetyl-CoA [nH=34]'] = Molecule.FromSmiles('CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')

        #mols['acetamide [nH=5]'] = Molecule.FromSmiles('CC(=O)N')
        #mols['acetamide [nH=4]'] = Molecule.FromSmiles('CC(=O)[NH-]')
        
        #mols['sparteine [nH=27]'] = Molecule.FromSmiles('[H][C@@]12CCCC[NH+]1C[C@@H]1C[C@H]2C[NH+]2CCCC[C@]12[H]')
        #mols['sparteine [nH=28]'] = Molecule.FromSmiles('[H][C@@]12CCCC[NH+]1C[C@@H]1C[C@H]2CN2CCCC[C@]12[H]')
        #mols['sparteine [nH=26]'] = Molecule.FromSmiles('[H][C@@]12CCCCN1C[C@@H]1C[C@H]2CN2CCCC[C@]12[H]')
        #mols['acetyl-CoA a'] = kegg.cid2mol(24)
        #mols['acetyl-CoA b'] = Molecule.FromSmiles("CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C")
        #mols['acetyl-CoA c'] = Molecule.FromSmiles("CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
        #mols['glycylglycine'] = Molecule.FromSmiles('C(C(=O)NCC(=O)[O-])[NH3+]')
        #mols['N-Acetylornithine'] = kegg.cid2mol(437)
        #mols['Sinapoyl-CoA'] = kegg.cid2mol(411)
        
        #smarts = 'C(=O)[N;H1;0]C'
        smarts = None
        
        for key, mol in mols.iteritems():
            #mol.SetTitle(key)
            print '-'*100
            print key
            if smarts:
                print mol.FindSmarts(smarts)
            try:
                decomposition = G.Mol2Decomposition(mol, ignore_protonations=False)
                print decomposition.ToTableString()
                print 'nH =', decomposition.Hydrogens()
                print 'z =', decomposition.NetCharge()
                print 'nMg = ', decomposition.Magnesiums()
                #for groupvec in decomposition.PseudoisomerVectors():
                #    print groupvec
                #    try:
                #        print G.groupvec2val(groupvec)
                #    except GroupMissingTrainDataError as e:
                #        print e.Explain(G)
            except GroupDecompositionError as e:
                print "Cannot decompose compound to groups: " + str(e)

            #try:
            #    pmap = G.GroupDecomposition2PseudoisomerMap(decomposition)
            #    print pmap
            #    dG0_tag = pmap.Transform(pH=7, pMg=14, I=0.1, T=298.15)
            #    print "dG0' = %.1f kJ/mol" % dG0_tag
            #except GroupMissingTrainDataError as e:
            #    print e.Explain(G)
            mol.Draw()
                
