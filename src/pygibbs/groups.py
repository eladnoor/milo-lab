#!/usr/bin/python

import csv
import logging
import os
import types
import json

from copy import deepcopy
import pylab
import sys
from pygibbs.thermodynamic_constants import R, default_pH, default_pMg, default_I, default_T, default_c0, dG0_f_Mg
from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException
from pygibbs.groups_data import Group, GroupsData
from pygibbs.group_vector import GroupVector
from pygibbs.pseudoisomers_data import PseudoisomersData, DissociationTable
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs import templates
from toolbox import util
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.sparse_kernel import SparseKernel, CplexNotInstalledError
from toolbox.molecule import Molecule

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
        
class GroupObservation(object):
    
    def __init__(self):
        self.groupvec = None
        self.dG0 = None
        self.name = None
        self.cid = None
        self.obs_type = None # can be 'formation', 'acid-base' or 'Mg'

    def NetCharge(self):
        return self.groupvec.NetCharge()

    def Hydrogens(self):
        return self.groupvec.Hydrogens()

    def Magnesiums(self):
        return self.groupvec.Magnesiums()

    def ToDatabase(self, db):
        db.Insert('group_observations', list(self.groupvec) +
                  [self.dG0, self.name, self.cid, self.obs_type])
        
class GroupObervationCollection(object):
    
    def __init__(self, db, html_writer, group_decomposer):
        self.kegg = Kegg.getInstance()
        self.db = db
        self.groups_data = group_decomposer.groups_data
        self.group_decomposer = group_decomposer
        self.observations = []
        self.n_groups = len(self.groups_data.all_groups)
        self.groupvec_titles = ["g%d" % i for i in xrange(self.n_groups)]
        self.cid2pmap_dict = {}
        self.cid_test_set = set()

        self.html_writer = html_writer
        if html_writer and html_writer.filename:
            self.FIG_DIR = os.path.basename(html_writer.filename).rsplit('.', 1)[0]
        else:
            self.FIG_DIR = ''
        
    def Add(self, groupvec, dG0, name, cid, obs_type):
        obs = GroupObservation()
        obs.groupvec = groupvec
        obs.dG0 = dG0
        obs.name = name
        obs.cid = cid
        obs.obs_type = obs_type
        self.observations.append(obs)

    def smiles2groupvec(self, smiles, id):
        mol = Molecule.FromSmiles(str(smiles))
        mol.RemoveHydrogens()
        mol.SetTitle(id)
        
        try:
            self.html_writer.write(mol.ToSVG())
        except (TypeError, AssertionError): # The Mg ions cannot be drawn by OASA 
            pass
        
        try:
            return self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
        except GroupDecompositionError as _e:
            logging.warning('Cannot decompose one of the compounds in the training set: %s, %s' % (id, smiles))
            self.html_writer.write('%s Could not be decomposed<br>\n' % smiles)

    def AddDissociationTable(self, cid, diss_table):
        name = self.kegg.cid2name(cid)
        logging.info("Reading pKa data for C%05d, %s:" % (cid, name))
        self.html_writer.write('<h3>C%05d - %s</h3>\n' % (cid, name))
        for key in diss_table:
            nH_above, nH_below, nMg_above, nMg_below = key
            ddG0, _ref = diss_table.ddGs[key]
            smiles_above, smiles_below = diss_table.smiles_dict[key]
            
            if not smiles_above or not smiles_below:
                continue
    
            if nH_above != nH_below:
                s = "nH = %d -> %d" % (nH_above, nH_below)
            elif nMg_above != nMg_below:
                s = "nMg = %d -> %d" % (nMg_above, nMg_below)
            logging.info("\t" + s)
            self.html_writer.write('%s: &#x394;G = %.2f<br>\n' 
                                   % (s, ddG0))
            self.html_writer.write('SMILES = %s >> %s<br>\n' % (smiles_below, smiles_above))
            decomposition_below = self.smiles2groupvec(smiles_below, 
                            "C%05d_b_H%d_Mg%d" % (cid, nH_below, nMg_below))
            decomposition_above = self.smiles2groupvec(smiles_above, 
                            "C%05d_a_H%d_Mg%d" % (cid, nH_above, nMg_above))
            if not decomposition_below or not decomposition_above:
                continue
            groupvec = decomposition_below.AsVector() - decomposition_above.AsVector()
            self.html_writer.write('<br>\nDecomposition = %s<br>\n' % str(groupvec))
            
            if nH_above != nH_below:
                if nH_above != nH_below-1:
                    raise Exception("a pKa must represent a difference of exactly 1 hydrogen")
                fullname = '%s [nH=%d]' % (name, nH_above)
                self.Add(groupvec, ddG0, fullname, cid, 'acid-base')
            elif nMg_above != nMg_below:
                if nMg_above != nMg_below-1:
                    raise Exception("a pK_Mg must represent a difference of exactly 1 magnesium")
                fullname = '%s [nMg=%d]' % (name, nMg_above)
                #self.Add(groupvec, ddG0, fullname, cid, 'Mg')
                # TODO: solve the problems with Mg species !!!
                # there seems to be a problem in decomposition of ADP:Mg (and maybe other forms that have Mg ions)

    def AddPseudoisomersData(self, ps_isomer):
        #name = str(ps_isomer)
        name = "%s (%d)" % (ps_isomer.name, ps_isomer.net_charge)
        logging.info('Verifying data for %s', name)
        self.html_writer.write("<h3>%s, %s</h3>\n" % (ps_isomer.name, ps_isomer.ref))

        if ps_isomer.dG0 == None:
            self.html_writer.write('No data for &#x394;G<sub>f</sub><br>\n')
            return

        if ps_isomer.Skip():
            self.html_writer.write('Compound marked as not to be used<br>\n')
            return
            
        self.html_writer.write('&#x394;G<sub>f</sub> = %.2f<br>\n' % ps_isomer.dG0)

        if ps_isomer.cid:
            self.AddPseudoisomer(ps_isomer.cid, ps_isomer.hydrogens, ps_isomer.net_charge,
                                ps_isomer.magnesiums, ps_isomer.dG0)

        if ps_isomer.Test():
            self.cid_test_set.add(ps_isomer.cid)
            self.html_writer.write('Compound marked to be used only for testing (not training)<br>\n')
            return
        
        elif ps_isomer.Train():
            self.html_writer.write('Compound marked to be used for training<br>\n')
        else:
            raise Exception('Unknown usage flag: %' % ps_isomer.use_for)

        if not ps_isomer.smiles:
            raise Exception("Cannot use compound '%s' for training if it lacks a SMILES string" % ps_isomer.name)
        try:
            self.html_writer.write('SMILES = %s<br>\n' % ps_isomer.smiles)
            mol = ps_isomer.MolNoH()
        except TypeError, e:
            logging.error(e)
            raise Exception('Invalid smiles: %s' % ps_isomer.smiles)

        mol.SetTitle(name)
        try:
            self.html_writer.write(mol.ToSVG())
        except (TypeError, IndexError, AssertionError):
            logging.warning('PyBel cannot draw the compound %s',  name)
            self.html_writer.write('WARNING: cannot draw this compound using PyBel\n')

        self.html_writer.write('<br>\n')

        try:
            decomposition = self.group_decomposer.Decompose(mol, strict=True)
        except GroupDecompositionError as e:
            logging.error('Cannot decompose one of the compounds in the training set: ' + str(mol))
            raise e
        
        groupvec = decomposition.AsVector()
        self.Add(groupvec, ps_isomer.dG0, name, ps_isomer.cid, 'formation')
        self.html_writer.write("Decomposition = %s<br>\n" % decomposition)
        
        gc_hydrogens, gc_charge = decomposition.Hydrogens(), decomposition.NetCharge()
        if ps_isomer.hydrogens != gc_hydrogens:
            s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, formula = %d' % (
                ps_isomer.hydrogens, gc_hydrogens)
            logging.error(s)
            self.html_writer.write(s + '<br>\n')
        if ps_isomer.net_charge != gc_charge:
            s = 'ERROR: Charge doesn\'t match: explicit = %d, formula = %d' % (
                ps_isomer.net_charge, gc_charge)
            logging.error(s)
            self.html_writer.write(s + '<br>\n')
        
    def ToDatabase(self):
        # This table is used only as an output for checking results
        # it is not used elsewhere in the code
        self.db.CreateTable('train_groups', 'name TEXT, nH INT, z INT, nMg INT')
        for group in self.groups_data.all_groups:
            self.db.Insert('train_groups', [group.name, group.hydrogens,
                                            group.charge, group.nMg])

        # the table 'group_observation' will contain all the observed data
        # that is used for training later
        titles = ['%s REAL' % t for t in self.groupvec_titles] + \
                 ['dG0 REAL', 'name TEXT', 'cid INT', 'obs_type TEXT']
        self.db.CreateTable('group_observations', ','.join(titles))

        self.db.CreateTable('pseudoisomer_observation', 
            'cid INT, name TEXT, protons INT, charge INT'
            ', nMg INT, dG0_f REAL, use_for TEXT')
        
        for obs in self.observations:
            obs.ToDatabase(self.db)
        self.db.Commit()
        
        for cid in self.cid2pmap_dict.keys():
            if cid in self.cid_test_set:
                use_for = 'test'
            else:
                use_for = 'train'
            for nH, z, nMg, dG0 in self.cid2pmap_dict[cid].ToMatrix():
                self.db.Insert('pseudoisomer_observation', 
                    [cid, self.kegg.cid2name(cid), nH, z, nMg, dG0, use_for])
    
    def SetPseudoisomerMap(self, cid, pmap):
        self.cid2pmap_dict[cid] = pmap

    def AddPseudoisomer(self, cid, nH, z, nMg, dG0):
        self.cid2pmap_dict.setdefault(cid, PseudoisomerMap())
        self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)

    def FromDatabase(self):
        self.observations = []
        for row in self.db.Execute("SELECT * FROM group_observations"):
            obs = GroupObservation()
            obs.groupvec = GroupVector(self.groups_data, row[0:self.n_groups])
            obs.dG0, obs.name, obs.cid, obs.obs_type = row[self.n_groups:]
            self.observations.append(obs)

        self.cid2pmap_dict = {}
        for row in self.db.Execute("SELECT * FROM pseudoisomer_observation"):
            cid, unused_name, nH, z, nMg, dG0, use_for = row
            self.AddPseudoisomer(cid, nH, z, nMg, dG0)
            if use_for == 'test':
                self.cid_test_set.add(cid)
        
    def GetRegressionData(self):
        """ 
            Returns the regression matrix and corresponding observation data 
            and names as a tuple: (A, b, names)
            Next step is to solve Ax = b
            
            Makes sure that there are no duplicates (two identical rows in A)
            by averaging the observed values.
        """
        
        hash2obs_vec = {}
        for obs in self.observations:
            h = str(obs.groupvec)
            hash2obs_vec.setdefault(h, []).append(obs)
        
        A = []
        b = []
        names = []
        for h, obs_vec in hash2obs_vec.iteritems():
            A.append(obs_vec[0].groupvec)
            
            all_cids = '; '.join(["%d" % obs.cid for obs in obs_vec])
            all_names = '; '.join([obs.name for obs in obs_vec])
            if obs_vec[0].obs_type == 'formation':
                names.append(all_names)
            else:
                names.append('pK:(%s) - [%s]' % (h, all_cids))
            b.append(pylab.mean([obs.dG0 for obs in obs_vec]))
            
        return pylab.matrix(A), pylab.array(b), names

class GroupContribution(Thermodynamics):    
    def __init__(self, db, html_writer=None):
        Thermodynamics.__init__(self)
        self.db = db

        if html_writer:
            self.html_writer = html_writer
        else:
            self.html_writer = NullHtmlWriter()

        self.kegg = Kegg.getInstance()
        self.bounds = deepcopy(self.kegg.cid2bounds)
            
        self.group_nullspace = None
        self.group_contributions = None
        self.obs_collection = None
            
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

    def read_training_data_formation(self, obs_fname="../data/thermodynamics/dG0.csv"):
        """
            Finds all the compounds which have a valid dG0 in the dG0.csv file,
            and adds the qualifying observations to obs_collection.
        """
        self.html_writer.write('<h2><a name=compounds>List of compounds for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.html_writer.write('Source File = %s<br>\n' % obs_fname)
        
        pdata = PseudoisomersData.FromFile(obs_fname)
        for ps_isomer in pdata:
            self.obs_collection.AddPseudoisomersData(ps_isomer)
        
        self.html_writer.write('</div>')
        
    def read_training_data_pKa(self):
        self.html_writer.write('<h2><a name=compounds>List of pKa for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)

        cid2DissociationTable = DissociationTable.ReadDissociationCsv()
        for cid, diss_table in cid2DissociationTable.iteritems():
            self.obs_collection.AddDissociationTable(cid, diss_table)
            
        self.html_writer.write('</div>')
            
    def load_training_data(self):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        self.obs_collection.FromDatabase()
        self.group_matrix, self.obs, self.mol_names = self.obs_collection.GetRegressionData()

    def train(self, FromFiles=True):
        self.obs_collection = GroupObervationCollection(self.db, 
            self.html_writer, self.group_decomposer)
        if FromFiles:
            self.read_training_data_pKa()
            self.read_training_data_formation()
            self.obs_collection.ToDatabase()
        else:
            self.obs_collection.FromDatabase()
            
        self.group_matrix, self.obs, self.mol_names = self.obs_collection.GetRegressionData()
        
        self.group_contributions, self.group_nullspace = self.RunLinearRegression()
        self.SaveContributionsToDB()
            
    def RunLinearRegression(self):
        group_contributions, _nullspace = LinearRegression.LeastSquares(self.group_matrix, self.obs)
        try:
            nullspace = SparseKernel(self.group_matrix).Solve()
        except CplexNotInstalledError:
            logging.warning("CPLEX is not installed on this system, using a non-sparse"
                            " method for describing the Kernel of the group matrix")
            nullspace = _nullspace
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
            self.db.Table2HTML(self.html_writer, 'train_group_matrix')
            self.db.Table2HTML(self.html_writer, 'train_observations')
            self.db.Table2HTML(self.html_writer, 'train_molecules')
            self.db.Table2HTML(self.html_writer, 'gc_contribution')
            
        self.html_writer.write('<table border="1">\n<tr>'
                        '<td>Name</td>'
                        '<td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td>'
                        '<td>Group Vector</td></tr>')
        for i in xrange(self.group_matrix.shape[0]):
            self.html_writer.write('<tr><td>%s</td><td>%.2f</td><td>%s</td></tr>' % \
                            (self.mol_names[i], self.obs[i], 
                             ' '.join(['%d' % self.group_matrix[i, j] for j in xrange(self.group_matrix.shape[1])])))
        self.html_writer.write('</table>')
        
        self.html_writer.write('<h2><a name="group_contrib">Group Contributions</a></h2>\n')
        self.html_writer.write('<table border="1">')
        self.html_writer.write('  <tr><td>#</td><td>Group Name</td><td>nH</td><td>charge</td><td>nMg</td><td>&#x394;<sub>gr</sub>G [kJ/mol]</td><td>&#x394;<sub>gr</sub>G\' [kJ/mol]</td><td>Appears in compounds</td></tr>\n')
        for i, group in enumerate(self.groups_data.all_groups):
            dG0_gr = self.group_contributions[i]
            dG0_gr_tag = dG0_gr + R*T*pylab.log(10)*pH*group.hydrogens
            compound_list_str = ' | '.join([self.mol_names[k] for k in pylab.find(self.group_matrix[:, i])])
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
            if self.mol_names[i][0:3] == 'pK:':
                continue
            logging.info('Cross validation, leaving-one-out: ' + self.mol_names[i])
            
            # leave out the row corresponding with observation 'i'
            subset = range(n_obs)
            subset.pop(i)
            
            group_contributions, nullspace = LinearRegression.LeastSquares(
                self.group_matrix[subset, :], self.obs[subset])
            
            if nullspace.shape[0] > self.group_nullspace.shape[0]:
                logging.warning('# example %d is not linearly dependent in the other examples' % i)
                deviations.append((None, self.mol_names[i], "%.1f" % self.obs[i], "-", "-", 
                                   "linearly independent example", "-"))
                continue
            
            orig_estimation = pylab.dot(self.group_matrix[i, :], self.group_contributions)
            estimation = pylab.dot(self.group_matrix[i, :], group_contributions)[0, 0]
            error = self.obs[i] - estimation
            
            val_names.append(self.mol_names[i])
            val_obs.append(self.obs[i])
            val_est.append(estimation)
            val_err.append(error)
            logging.info('Error = %.1f' % error)
            deviations.append((abs(error), self.mol_names[i], 
                               "%.1f" % self.obs[i], 
                               "%.1f" % estimation,
                               "%.1f" % error, "",
                               "%.1f" % orig_estimation))
        
        logging.info("writing the table of estimation errors for each compound")
        self.html_writer.write('<h2><a name=compounds>Cross Validation Table</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)

        self.html_writer.write('<div><b>std(error) = %.2f kJ/mol</b></div>\n' % pylab.std(val_err))
    
        self.html_writer.write('<div><b>rmse(pred) = %.2f kJ/mol</b></div>\n' % util.calc_rmse(val_obs, val_est))

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

    def cid2SourceString(self, cid):
        """ override the method in Thermodynamics """
        return "Group Contribution"
        
    def cid2PseudoisomerMap(self, cid):
        """
            Returns a PseudoisomerMap according to a given CID.
        """
        # if the entire KEGG compound database has been computed in advance, use the cached data
        if cid == 80: # H+
            pmap = PseudoisomerMap()
            pmap.Add(0, 0, 0, 0)
            return pmap

        try:
            mol = self.kegg.cid2mol(cid)
            mol.title = "C%05d" % cid
            return self.Mol2PseudoisomerMap(mol, ignore_protonations=True)
        except KeggParseException as e:
            error = "Cannot determine molecular structure: " + str(e)
        except GroupDecompositionError as e:
            error = "Cannot decompose into groups: " + str(e)
        except GroupMissingTrainDataError as e:
            error = "Cannot calculate PseudoisomerMap: " + str(e)

        #if cid in self.obs_collection.cid2pmap_dict:
        #    return self.obs_collection.cid2pmap_dict[cid]
        raise MissingCompoundFormationEnergy(error, cid)

    def estimate_pKa_keggcid(self, cid, charge, T=default_T):
        """
            Estimates the pKa of the compound.
            Ka is the equilibrium constant for the protonation reaction
            from the pseudoisomer with 'charge' to the pseudoisomer with 'charge'+1 
        """
        dG0_p0 = None
        dG0_p1 = None
        for dG0, unused_nH, z, unused_nMg in self.cid2PseudoisomerMap(cid).ToMatrix():
            if z == charge:
                dG0_p0 = dG0
            elif z == charge + 1:
                dG0_p1 = dG0
        
        if dG0_p0 == None:
            raise GroupMissingTrainDataError(
                "cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge))
        if dG0_p1 == None:
            raise GroupMissingTrainDataError(
                "cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge + 1))
        
        return (dG0_p0 - dG0_p1) / (R * T * pylab.log(10))

    def estimate_dG0(self, mol, pH=default_pH, pMg=default_pMg,
                     I=default_I, T=default_T):
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).
            
            T in Kelvin
        """
        pmap = self.Mol2PseudoisomerMap(mol)
        return pmap.TransformMatrix(pH, pMg, I, T)
            
    def estimate_dG0_keggcid(self, cid, pH=default_pH, pMg=default_pMg,
                             I=default_I, T=default_T, most_abundant=False): # T = temperature in K
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).            
        """
        pmap = self.cid2PseudoisomerMap(cid)
        return pmap.TransformMatrix(pH, pMg, I, T)
    
    def estimate_dG0_reaction(self, sparse_reaction, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, most_abundant=False):
        """
            Calculate the dG0 (standard Gibbs free energy change)
            pH and I can be either floats or lists of values. If both are floats, returns a float.
            If either pH and/or I are lists, returns a matrix.
            T must be a float.
            
            dG = dG0 + R T \sum_i{s_i * ln([C_i])}
            
            s_i - stoichiometric coeff of compound i
            [C_i] - concentration of compound C_i
        """
        
        pmaps = [self.cid2PseudoisomerMap(cid) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()

        if type(pH) != types.ListType and type(I) != types.ListType:
            dG0_vector = [pmap.Transform(pH, pMg, I, T, most_abundant) for pmap in pmaps]
            dG0 = pylab.dot(stoichiometry_vector, dG0_vector)
            return dG0
        else:
            if type(pH) != types.ListType:
                pH = [float(pH)]
            if type(I) != types.ListType:
                I = [float(I)]
            
            dG0_matrix = pylab.zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_vector = [pmap.Transform(pH[i], pMg, I[j], T, most_abundant) for pmap in pmaps]
                    dG0_matrix[i, j] = pylab.dot(stoichiometry_vector, dG0_vector)
            
            return dG0_matrix

    def estimate_dG_reaction(self, sparse_reaction, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, c0=default_c0, media=None, most_abundant=False):
        """
            standard = False means to use the known concentrations of reactants as well
        """
        
        dG0 = self.estimate_dG0_reaction(sparse_reaction, pH, pMg, I, T, most_abundant)
        concentration_vector = [pylab.log(self.get_concentration(cid, c0, media)) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()
        
        return dG0 + R * T * pylab.dot(stoichiometry_vector, concentration_vector)        

    def estimate_dG_keggrid(self, rid, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, c0=default_c0, media=None, most_abundant=False):
        """
            Returns the transformed Gibbs free energy change of a reaction according to its RID.
            Can set the pH, I and T to non-standard values.
            When media == None, it means we should use standard conditions (i.e. dG0).
            When media == 'glucose' (for example), it uses the concentrations measured for growth on glucose media.
        """
        sparse_reaction = self.kegg.rid2sparse_reaction(rid) 
        try:
            return self.estimate_dG_reaction(sparse_reaction, pH, pMg, I, T, c0, media, most_abundant)
        except KeyError as e:
            raise KeyError("R%05d contains a compound which cannot be used\n" % rid + str(e))

    def estimate_dG0_reaction_formula(self, formula, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, most_abundant=False):
        sparse_reaction = self.kegg.formula_to_sparse(formula)
        return self.estimate_dG0_reaction(sparse_reaction, pH, pMg, I, T, most_abundant)

    def estimate_dG_reaction_formula(self, formula, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, media=None, most_abundant=False):
        sparse_reaction = self.kegg.formula_to_sparse(formula)
        return self.estimate_dG_reaction(sparse_reaction, pH, pMg, I, T, media, most_abundant)
        
    def cid2groupvec(self, cid):
        try:
            return self.get_groupvec(self.kegg.cid2mol(cid))
        except GroupDecompositionError:
            raise GroupDecompositionError("Unable to decompose %s (C%05d) into groups" % (self.kegg.cid2name(cid), cid))
            
    def rid2groupvec(self, rid):
        sparse_reaction = self.kegg.rid2sparse_reaction(rid) 
        group_matrix = pylab.matrix([self.cid2groupvec(cid) for cid in sparse_reaction.keys()])
        stoichiometry_vector = pylab.matrix(sparse_reaction.values())
        total_groupvec = pylab.dot(stoichiometry_vector, group_matrix)
        return total_groupvec.tolist()[0]

    def write_cid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(['cid'] + self.all_group_names)
        for cid in self.kegg.get_all_cids_with_inchi():
            try:
                groupvec = self.cid2groupvec(cid)
                csv_file.writerow([cid] + groupvec)
            except GroupDecompositionError as e:
                print str(e)
                continue
            except KeggParseException as e:
                print str(e)
                continue
            
    def write_rid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(['rid'] + self.all_group_names)
        group_names = self.all_group_names
        
        groupstr_to_counter = {}
        for rid in self.kegg.get_all_rids():
            try:
                groupvec = self.rid2groupvec(rid)
                csv_file.writerow([rid] + groupvec)
                groupstr = ""
                for i in range(len(group_names)):
                    if (groupvec[i] != 0):
                        groupstr += "%s : %d, " % (group_names[i], groupvec[i])
                groupstr_to_counter[groupstr] = groupstr_to_counter.get(groupstr, 0) + 1
            except GroupDecompositionError as e:
                print "R%05d: Cannot decompose at least one of the compounds\n" % rid + str(e)
            except KeggParseException as e:
                print str(e)
            except KeyError as e:
                print "R%05d: Cannot locate the reaction or one of the compounds in kegg\n" % rid + str(e)
        
        gstr_hist = [(counter, groupstr) for (groupstr, counter) in groupstr_to_counter.iteritems()]
        gstr_hist.sort(reverse=True)
        f = open("../res/groupstr_hist.txt", "w")
        for (counter, groupstr) in gstr_hist:
            f.write("%5d ; " % counter + groupstr + "\n")
        f.close()
        
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
                self.db.Insert('compound_abundance', [cid, "glucose", float(row['Glycerol'])])
            except ValueError:
                pass
            try:
                self.db.Insert('compound_abundance', [cid, "glucose", float(row['Acetate'])])
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

#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
    
if __name__ == '__main__':
    db = SqliteDatabase('../res/gibbs.sqlite')
    if len(sys.argv) < 2:
        html_writer = HtmlWriter('../res/groups.html')
        G = GroupContribution(db=db, html_writer=html_writer)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.train()
        G.write_regression_report()
        G.analyze_training_set()
        logging.info("Estimating formation energies for all KEGG. Please be patient for a few minutes...")
        G.ToDatabase(db, table_name='gc_pseudoisomers', error_table_name='gc_errors')
    else:
        G = GroupContribution(db=db)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.init()

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
                decomposition = G.Mol2Decomposition(mol, ignore_protonations=True)
                print decomposition.ToTableString()
                #for groupvec in decomposition.PseudoisomerVectors():
                #    print groupvec
                #    try:
                #        print G.groupvec2val(groupvec)
                #    except GroupMissingTrainDataError as e:
                #        print e.Explain(G)
                pmap = G.GroupDecomposition2PseudoisomerMap(decomposition)
                print pmap
                dG0 = pmap.Transform()
                print dG0
            except GroupDecompositionError as e:
                print "Cannot decompose compound to groups: " + str(e)
            except GroupMissingTrainDataError as e:
                print e.Explain(G)
            mol.Draw()
                
