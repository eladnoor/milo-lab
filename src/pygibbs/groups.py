#!/usr/bin/python

import csv
import logging
import openbabel
import os
import pybel
import types
import json

from copy import deepcopy
import pylab
from pygibbs.thermodynamic_constants import R, default_pH, default_pMg, default_I, default_T, default_c0, dG0_f_Mg
from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.kegg import KeggParseException, Kegg
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.groups_data import Group, GroupsData
from pygibbs import templates
from toolbox import util
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from pygibbs.pseudoisomers_data import PseudoisomersData
from pygibbs.pseudoisomer import PseudoisomerMap
import sys

class GroupContributionError(Exception):
    pass
    
class GroupMissingTrainDataError(Exception):
    def __init__(self, value, kernel_rows=[]):
        self.value = value
        self.kernel_rows = kernel_rows
    def __str__(self):
        if (type(self.value) == types.StringType):
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
        Thermodynamics.__init__(self)
        self.db = db

        if html_writer:
            self.FIG_DIR = os.path.basename(html_writer.filename).rsplit('.', 1)[0]
            self.html_writer = html_writer
        else:
            self.FIG_DIR = ''
            self.html_writer = NullHtmlWriter()

        if kegg:
            self._kegg = kegg
        else:
            self._kegg = Kegg(self.db)
        self.bounds = deepcopy(self._kegg.cid2bounds)
        self.hatzi = Hatzi()
        self.hatzi.bounds = self.bounds
            
        self.override_gc_with_measurements = True
        self.cid2pmap_dict = None
        self.group_nullspace = None
        self.group_contributions = None
            
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
        self.load_cid2pmap()
        
    def quick_init(self, cid_list=None, html_fname=None):
        self.load_groups()
        self.LoadContributionsFromDB()
        self.load_concentrations()
        self.load_training_data()
        self.quick_calculate_pmaps(cid_list, html_fname)

    def quick_calculate_pmaps(self, cid_list=None, html_fname=None):
        """
            Use this method only if you want a quick recalculation of some
            of the CIDs. Note that it does not generate any logs or reports,
            or write the results to the database.
            Please use 'save_cid2pmap' in most cases, and then load the pmaps
            from the database using 'load_cid2pmap'.
        """
        data = {}
        data['compounds'] = []
        self.cid2pmap_dict = {}
        
        cid_list = cid_list or self.kegg().get_all_cids()
        logging.info('Recalculating formation energies for %d compounds' % len(cid_list))
        for cid in cid_list:
            comp = self.kegg().cid2compound(cid)

            cdict = {}
            cdict['cid'] = cid
            cdict['compound'] = comp
            cdict['measured_pmap'] = self.cid2pmap_obs.get(cid, None)
            cdict['estimated_pmap'] = None
            
            if comp.inchi:
                try:
                    decomposition = self.Mol2Decomposition(comp.get_mol(),
                                                           ignore_protonations=True)
                    cdict['decomposition'] = decomposition
                    pmap = self.GroupDecomposition2PseudoisomerMap(decomposition)
                    cdict['estimated_pmap'] = pmap
                except (KeggParseException, GroupDecompositionError, 
                        GroupMissingTrainDataError):
                    pass
            
            pmap = cdict['measured_pmap'] or cdict['estimated_pmap']
            if pmap:
                self.cid2pmap_dict[cid] = pmap
            data['compounds'].append(cdict)

        if html_fname:
            templates.render_to_file('kegg_pmaps.html', data, html_fname)
                
    def save_cid2pmap(self):
        self.cid2pmap_dict = {}
        logging.info('Calculating the table of chemical formation energies for all KEGG compounds.')
        self.db.CreateTable('gc_cid2prm', 'cid INT, nH INT, z INT, nMg INT, dG0 REAL, estimated BOOL')
        self.db.CreateTable('gc_cid2error', 'cid INT, error TEXT')

        compounds = []

        for cid in self.kegg().get_all_cids():
            cdict = {'cid': cid, 'measured_pmap': None,
                     'estimated_pmap': None, 'compound': None}
            
            if cid % 100 == 1:
                logging.info('Saving KEGG Compound C%05d' % cid)
            
            # If the compound is measured:
            if (cid in self.cid2pmap_obs):
                pmap = self.cid2pmap_obs[cid]
                cdict['measured_pmap'] = pmap
                for (nH, z, nMg, dG0) in pmap.ToMatrix():
                    self.db.Insert('gc_cid2prm', [cid, nH, z, nMg, dG0, False])

            # Try to also estimate the dG0_f using Group Contribution:
            comp = self.kegg().cid2compound(cid)
            cdict['compound'] = comp
            error_str = None
            if not comp.inchi:
                error_str = 'no InChI exists'
            else:
                try:
                    decomposition = self.Mol2Decomposition(comp.get_mol(),
                                                           ignore_protonations=True)
                    cdict['decomposition'] = decomposition
                    pmap = self.GroupDecomposition2PseudoisomerMap(decomposition)
                    cdict['estimated_pmap'] = pmap
                    for (nH, z, nMg, dG0) in pmap.ToMatrix():
                        self.db.Insert('gc_cid2prm', [cid, int(nH), int(z), int(nMg), dG0, True])
                except KeggParseException:
                    error_str = 'cannot determine molecular structure'
                except GroupDecompositionError:
                    error_str = 'cannot decompose into groups'
                except GroupMissingTrainDataError as e:
                    error_str = e.Explain(self)
            compounds.append(cdict)
            self.db.Insert('gc_cid2error', [cid, error_str])
        
        templates.render_to_file('kegg_pmaps.html', {'compounds': compounds},
                                 '../res/kegg_pmaps.html')
        logging.info('Writing the formation energies to the database')
        self.db.Commit()
        logging.info('DONE!')

    def load_cid2pmap(self):
        self.cid2pmap_dict = {}
        self.cid2source_string = {}

        # Now load the data into the cid2pmap_dict:
        for row in self.db.Execute("SELECT cid, nH, z, nMg, dG0 from gc_cid2prm WHERE estimated == 1;"):
            cid, nH, z, nMg, dG0 = row
            self.cid2pmap_dict.setdefault(cid, PseudoisomerMap())
            self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)
            self.cid2source_string[cid] = 'Group Contribution'

        cid2pmap_obs = {} # observed formation energies
        for row in self.db.Execute("SELECT cid, nH, z, nMg, dG0 from gc_cid2prm WHERE estimated == 0;"):
            cid, nH, z, nMg, dG0 = row
            cid2pmap_obs.setdefault(cid, PseudoisomerMap())
            cid2pmap_obs[cid].Add(nH, z, nMg, dG0)

        # add the observed data to the cid2pmap_dict (and override the estimated data if required)
        for cid in cid2pmap_obs.keys():
            if (self.override_gc_with_measurements or cid not in self.cid2pmap_dict):
                self.cid2pmap_dict[cid] = cid2pmap_obs[cid]
                self.cid2source_string[cid] = 'Observed'
        
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
            h['source'] = self.cid2source_string.get(row['cid'], None)
            h['error'] = row['error']
            h['species'] = []
            try:
                for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(row['cid']).ToMatrix():
                    h['species'].append({"nH":nH, "z":z, "nMg":nMg, "dG0_f":dG0})
            except MissingCompoundFormationEnergy:
                pass

            formations.append(h)

        json_file = open(json_fname, 'w')
        json_file.write(json.dumps(formations, indent=4))
        json_file.close()
        
    def train(self, obs_fname, use_dG0_format=False):
        l_groupvec_pKa, l_deltaG_pKa, l_names_pKa = self.read_training_data_pKa()
        #l_groupvec_pKa, l_deltaG_pKa, l_names_pKa = [], [], []

        if use_dG0_format:
            l_groupvec_formation, l_deltaG_formation, l_names_formation = self.read_training_data_dG0(obs_fname)
        else:
            l_groupvec_formation, l_deltaG_formation, l_names_formation = self.read_training_data(obs_fname)
        
        l_groupvec = l_groupvec_pKa + l_groupvec_formation
        l_deltaG = l_deltaG_pKa + l_deltaG_formation
        l_names = l_names_pKa + l_names_formation
        self.save_training_data(l_groupvec, l_deltaG, l_names)

        self.group_contributions, self.group_nullspace = \
            self.linear_regression_train()
        self.SaveContributionsToDB()
            
    def linear_regression_train(self):
        self.load_training_data()
        
        #group_contributions, nullspace = LinearRegression.LeastSquares(
        #    self.group_matrix, self.obs)
        group_contributions, nullspace = LinearRegression.SolveLinearSystem(
            self.group_matrix, self.obs)
        return list(group_contributions.flat), nullspace
    
    def load_groups(self, group_fname=None):
        if group_fname:
            self.groups_data = GroupsData.FromGroupsFile(group_fname)
            self.groups_data.ToDatabase(self.db)
        else:
            self.groups_data = GroupsData.FromDatabase(self.db)
        self.group_decomposer = GroupDecomposer(self.groups_data)
        self.dissociation = DissociationConstants(self.db, self.html_writer,
                                                  self.kegg(), self.group_decomposer)

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
        if (cid == 1): # the concentration of water must always be 1
            return 1
        if (media == None):
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
        
    def cid2pseudoisomers(self, cid):
        try:
            comp = self.kegg().cid2compound(cid)
            return self.get_pseudoisomers(comp.get_mol())
        except GroupDecompositionError:
            return [(self.kegg().cid2num_hydrogens(cid), self.kegg().cid2charge(cid), 0)]
        except KeggParseException:
            return [(0, 0, 0)]
    
    def mol2inchi(self, mol):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "inchi")
        return obConversion.WriteString(mol.OBMol).strip()        
    
    def read_training_data_dG0(self, obs_fname):
        """
            Finds all the compounds which have a valid dG0 in the dG0.csv file,
            and generates a regression matrix using the KEGG groups for these compounds.
            Return values is a tuple (X, dG_obs, names) where:
            X           - is the group regression matrix
            dG_obs      - is the observed dG0 vector.
            names       - is the list of compound names.
        """
        l_groupvec = []
        l_deltaG = []
        l_names = [] # is the list of Molecules (pybel class) used for the regression
        self.inchi2val_obs = {}
        self.cid2pmap_obs = {}
        self.cid_test_set = set()
        
        self.html_writer.write('<h2><a name=compounds>List of compounds for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.html_writer.write('Source File = %s<br>\n' % obs_fname)
        
        pdata = PseudoisomersData.FromFile(obs_fname)
        
        counter = 0
        for ps_isomer in pdata:
            name = str(ps_isomer)
            logging.info('Verifying data for %s', name)
            
            self.html_writer.write("<h3>%s, %s</h3>\n" % (ps_isomer.name, ps_isomer.ref))

            if ps_isomer.dG0 == None:
                self.html_writer.write('No data for &#x394;G<sub>f</sub><br>\n')
                continue

            if ps_isomer.Skip():
                self.html_writer.write('Compound marked as not to be used<br>\n')
                continue
                
            self.html_writer.write('&#x394;G<sub>f</sub> = %.2f<br>\n' % ps_isomer.dG0)

            if ps_isomer.cid:
                self.cid2pmap_obs.setdefault(ps_isomer.cid, PseudoisomerMap())
                self.cid2pmap_obs[ps_isomer.cid].Add(ps_isomer.hydrogens, ps_isomer.net_charge,
                                                     ps_isomer.magnesiums, ps_isomer.dG0)

            if ps_isomer.Test():
                self.cid_test_set.add(ps_isomer.cid)
                self.html_writer.write('Compound marked to be used only for testing (not training)<br>\n')
                continue
            
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

            inchi = self.mol2inchi(mol)
            self.html_writer.write('INCHI = %s<br>\n' % inchi)
            self.inchi2val_obs[inchi] = ps_isomer.dG0
            mol.title = name
            
            img_fname = self.FIG_DIR + '/train_%05d.png' % counter
            counter += 1
            try:
                self.html_writer.embed_molecule_as_png(mol, img_fname)
            except (TypeError, IndexError, AssertionError):
                logging.warning('PyBel cannot draw the compound %s',  name)
                self.html_writer.write('WARNING: cannot draw this compound using PyBel\n')

            self.html_writer.write('<br>\n')
    
            try:
                decomposition = self.group_decomposer.Decompose(mol, strict=True)
            except GroupDecompositionError as e:
                logging.error('Cannot decompose one of the compounds in the training set: ' + mol.title)
                raise e
            
            groupvec = decomposition.AsVector()
            l_names.append(name)
            l_groupvec.append(groupvec)
            l_deltaG.append(ps_isomer.dG0)
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
                
        self.html_writer.write('</div>')
        return (l_groupvec, l_deltaG, l_names)        
        
    def read_training_data(self, obs_fname):
        l_groupvec = []
        l_deltaG = []
        l_names = [] # is the list of Molecules (pybel class) used for the regression
        self.inchi2val_obs = {}
        self.cid2pmap_obs = {}
        
        self.html_writer.write('<h2><a name=compounds>List of compounds for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.html_writer.write('Source File = %s<br>\n' % obs_fname)
        Km_csv = csv.reader(open(obs_fname, 'r'))
        Km_csv.next()
        
        for row in Km_csv:
            (cid, obs_val) = row
            if (cid == "" or obs_val == ""):
                continue
            cid = int(cid)
            obs_val = float(obs_val)
            self.cid2pmap_obs[cid].Add(0, 0, 0, obs_val)
            try:
                mol = self.kegg().cid2mol(cid)
                self.html_writer.write('<h3><a href="%s">C%05d - %s</a></h3>\n' % (self.kegg().cid2link(cid), cid, self.kegg().cid2name(cid)))
                self.html_writer.write('Observed value = %f <br>\n' % obs_val)
                
                inchi = self.mol2inchi(mol)
                self.html_writer.write('INCHI = %s<br>\n' % inchi)
                self.inchi2val_obs[inchi] = obs_val
                
                groupvec = self.group_decomposer.Decompose(mol).AsVector()
                l_names.append(mol)
                l_groupvec.append(groupvec)
                l_deltaG.append(obs_val)
                self.html_writer.write('Decomposition = %s <br>\n' % self.get_decomposition_str(mol))
            except GroupDecompositionError:
                self.html_writer.write('Could not be decomposed<br>\n')
            except KeyError:
                self.html_writer.write('Compound has no INCHI in KEGG<br>\n')

        self.html_writer.write('</div>')
        return (l_groupvec, l_deltaG, l_names)        

    def read_training_data_pKa(self):
        
        def smiles2groupvec(smiles, id):
            mol = pybel.readstring('smiles', str(smiles))
            mol.removeh()
            mol.title = id
            self.html_writer.embed_molecule_as_png(mol, 'dissociation_constants/%s.png' % id)
            try:
                return self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
            except GroupDecompositionError as _e:
                logging.warning('Cannot decompose one of the compounds in the training set: %s, %s' % (id, smiles))
                self.html_writer.write('%s Could not be decomposed<br>\n' % smiles)
                #raise e
        
        l_groupvec = []
        l_deltaG = []
        l_names = [] # is the list of Molecules' names

        self.dissociation.LoadValuesToDB()
        self.html_writer.write('<h2><a name=compounds>List of pKa for training</a>')
        div_id = self.html_writer.insert_toggle()
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        for cid, T, nH_below, nH_above, smiles_below, smiles_above, pKa in self.db.Execute(
                "SELECT cid, T, nH_below, nH_above, smiles_below, smiles_above, pKa FROM pKa"):
            if not smiles_above or not smiles_below or not cid:
                continue
            
            logging.info("Reading pKa data for C%05d, %s, nH=%d -> nH=%d" % \
                         (cid, self.kegg().cid2name(cid), nH_below, nH_above))
            self.html_writer.write('<h3>C%05d - %s - nH=%d -> nH=%d</h3>\n' % \
                                   (cid, self.kegg().cid2name(cid), nH_below, nH_above))
            
            dG0 = R*T*pylab.log(10)*pKa
            self.html_writer.write('pKa = %.2f, T = %.2f \n</br>' % (pKa, T))
            self.html_writer.write('&#x394;G<sub>p</sub> = %.2f<br>\n' % dG0)
            
            self.html_writer.write('SMILES = %s >> %s<br>\n' % (smiles_below, smiles_above))
            decomposition_below = smiles2groupvec(smiles_below, "C%05d_b_H%d" % (cid, nH_below))
            decomposition_above = smiles2groupvec(smiles_above, "C%05d_a_H%d" % (cid, nH_above))
            if not decomposition_below or not decomposition_above:
                continue
            groupvec = decomposition_above.AsVector() - decomposition_below.AsVector()
            self.html_writer.write('<br>\nDecomposition = %s<br>\n' % str(groupvec))

            try:
                i = l_groupvec.index(groupvec)
                l_deltaG[i].append(dG0)
                l_names[i] += ', C%05d [%d->%d]' % (cid, nH_above, nH_below)
            except ValueError:
                l_groupvec.append(groupvec)
                l_deltaG.append([dG0])
                l_names.append('pKa (%s): C%05d [%d->%d]' % \
                               (str(groupvec), cid, nH_above, nH_below))
        
        for i in xrange(len(l_groupvec)):
            l_deltaG[i] = pylab.mean(l_deltaG[i])
        
        self.html_writer.write('</div>')
        return (l_groupvec, l_deltaG, l_names)        

    def save_training_data(self, l_groupvec, l_deltaG, l_names):
        n_obs = len(l_deltaG)
        if not n_obs:
            raise Exception("No observations have been given for training")
        
        n_groups = len(l_groupvec[0])
        
        self.db.CreateTable('train_group_matrix', ",".join(["g%d REAL" % i for i in range(n_groups)]))
        for j in range(n_obs):
            self.db.Insert('train_group_matrix', l_groupvec[j])
        
        self.db.CreateTable('train_observations', 'obs REAL')
        for i in range(n_obs):
            self.db.Insert('train_observations', [l_deltaG[i]])
        
        self.db.CreateTable('train_groups', 'name TEXT, nH INT, z INT, nMg INT')
        for group in self.groups_data.all_groups:
            self.db.Insert('train_groups', [group.name, group.hydrogens,
                                            group.charge, group.nMg])

        self.db.CreateTable('train_molecules', 'name TEXT')
        for name in l_names:
            self.db.Insert('train_molecules', [name])

        self.db.Commit()
            
    def load_training_data(self):
        X = []
        for row in self.db.Execute("SELECT * FROM train_group_matrix"):
            X.append(list(row))
        self.group_matrix = pylab.array(X)

        y = []
        for row in self.db.Execute("SELECT * FROM train_observations"):
            y.append(row[0])
        self.obs = pylab.array(y)
        
        self.mol_names = []
        for row in self.db.Execute("SELECT name FROM train_molecules"):
            self.mol_names.append(row[0])

    def groupvec2val(self, groupvec):
        if self.group_contributions == None or self.group_nullspace == None:
            raise Exception("You need to first Train the system before using it to estimate values")

        v = abs(pylab.dot(self.group_nullspace, pylab.array(groupvec)))
        k_list = [i for i in pylab.find(v > 1e-10)]
        if k_list:
            raise GroupMissingTrainDataError("can't estimate because the input "
                "is not orthogonal to the kernel", k_list)
        return pylab.dot(groupvec, self.group_contributions)
    
    def estimate_val(self, mol):
        try:
            groupvec = self.get_groupvec(mol)
        except GroupDecompositionError as e:
            inchi = self.mol2inchi(mol)
            if (inchi in self.inchi2val_obs):
                return self.inchi2val_obs[inchi]
            else:
                raise e
        
        return self.groupvec2val(groupvec)

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
                            (self.mol_names[i], self.obs[i], ' '.join(['%d' % x for x in self.group_matrix[i, :]])))
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
        self.html_writer.write('</h2><div id="%s" style="display:none">\n' % div_id)
        self.db.Table2HTML(self.html_writer, 'gc_nullspace')
        self.html_writer.write('</div>\n')

    def analyze_training_set(self):
        n_obs = self.group_matrix.shape[0]
        orig_est = []
        val_names = []
        val_obs = []
        val_est = []
        val_err = []
        deviations = []
        
        for i in xrange(n_obs):
            # skip the cross-validation of the pKa values since group
            # contribution is not meant to give any real prediction for pKas
            # except the mean of the class of pKas.
            if self.mol_names[i][0:3] == 'pKa':
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
            estimation = pylab.dot(self.group_matrix[i, :], group_contributions)[0]
            error = self.obs[i] - estimation
            
            val_names.append(self.mol_names[i])
            orig_est.append(orig_estimation)
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

    def kegg(self):
        return self._kegg        

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

    def Mol2Decomposition(self, mol, ignore_protonations=False):
        try:
            return self.group_decomposer.Decompose(
                mol, ignore_protonations, strict=True)
        except GroupDecompositionError as e:
            raise GroupDecompositionError(str(e) + "\n" + mol.title + "\n")

    def Mol2PseudoisomerMap(self, mol, ignore_protonations=False):
        decomposition = self.Mol2Decomposition(mol, ignore_protonations)
        return self.GroupDecomposition2PseudoisomerMap(decomposition)
    
    def GroupDecomposition2PseudoisomerMap(self, decomposition):
        all_groupvecs = decomposition.PseudoisomerVectors()
        if not all_groupvecs:
            raise GroupDecompositionError('Found no pseudoisomers for %s'
                                          % mol.title)

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
                                             decomposition.mol.title, kernel_rows)            

        pmap.Squeeze()
        return pmap

    def cid2PseudoisomerMap(self, cid, use_cache=True):
        """
            returns a list of 3-tuples of dG0 (untransformed), nH and z.
            Each tuple represents one of the pseudoisomers.
        """
        # if the entire KEGG compound database has been computed in advance, use the cached data
        if use_cache and self.cid2pmap_dict != None:
            if cid in self.cid2pmap_dict:
                return self.cid2pmap_dict[cid]
            else:
                raise MissingCompoundFormationEnergy("Formation energy cannot be determined using group contribution", cid)

        if (cid == 80): # H+
            pmap = PseudoisomerMap()
            pmap.Add(0, 0, 0, 0)
            return pmap

        if (cid in self.cid2pmap_obs and (self.override_gc_with_measurements or cid in self.cid_test_set)):
            pmap = PseudoisomerMap()
            pmap.AddAll(self.cid2pmap_obs[cid])
            return pmap
        else:
            try:
                mol = self.kegg().cid2mol(cid)
            except KeggParseException as e:
                raise MissingCompoundFormationEnergy("Cannot determine molecular structure: " + str(e), cid)
            
            mol.title = "C%05d" % cid
            try:
                return self.Mol2PseudoisomerMap(mol)
            except GroupContributionError as e:
                if cid in self.cid2pmap_obs:
                    return self.cid2pmap_obs[cid]
                else:
                    raise MissingCompoundFormationEnergy(str(e), cid)
    
    def estimate_pKa_keggcid(self, cid, charge, T=default_T):
        """
            Estimates the pKa of the compound.
            Ka is the equilibrium constant for the protonation reaction
            from the pseudoisomer with 'charge' to the pseudoisomer with 'charge'+1 
        """
        dG0_p0 = None
        dG0_p1 = None
        for (dG0, unused_nH, z, unused_nMg) in self.cid2PseudoisomerMap(cid).ToMatrix():
            if (z == charge):
                dG0_p0 = dG0
            elif (z == charge + 1):
                dG0_p1 = dG0
        
        if (dG0_p0 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge))
        if (dG0_p1 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge + 1))
        
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
        sparse_reaction = self.kegg().rid2sparse_reaction(rid) 
        try:
            return self.estimate_dG_reaction(sparse_reaction, pH, pMg, I, T, c0, media, most_abundant)
        except KeyError as e:
            raise KeyError("R%05d contains a compound which cannot be used\n" % rid + str(e))

    def estimate_dG0_reaction_formula(self, formula, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, most_abundant=False):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.estimate_dG0_reaction(sparse_reaction, pH, pMg, I, T, most_abundant)

    def estimate_dG_reaction_formula(self, formula, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, media=None, most_abundant=False):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.estimate_dG_reaction(sparse_reaction, pH, pMg, I, T, media, most_abundant)
        
    def cid2groupvec(self, cid):
        try:
            return self.get_groupvec(self.kegg().cid2mol(cid))
        except GroupDecompositionError:
            raise GroupDecompositionError("Unable to decompose %s (C%05d) into groups" % (self.kegg().cid2name(cid), cid))
            
    def rid2groupvec(self, rid):
        sparse_reaction = self.kegg().rid2sparse_reaction(rid) 
        group_matrix = pylab.matrix([self.cid2groupvec(cid) for cid in sparse_reaction.keys()])
        stoichiometry_vector = pylab.matrix(sparse_reaction.values())
        total_groupvec = pylab.dot(stoichiometry_vector, group_matrix)
        return total_groupvec.tolist()[0]

    def analyze_all_kegg_compounds(self, pH=[default_pH], pMg=default_pMg, I=[default_I], T=default_T, most_abundant=False):
        self.db.CreateTable('dG0_f', 'cid INT, pH REAL, pMg REAL, I REAL, T REAL, dG0 REAL')
        self.db.CreateIndex('dG0_f_idx', 'dG0_f', 'cid, pH, I, T', unique=True)

        self.html_writer.write('<h2><a name=kegg_compounds>&#x394;G<sub>f</sub> of KEGG compounds:</a></h2>')
        for cid in self.cid2pmap_dict.keys():
            self.html_writer.write('<p>\n')
            self.html_writer.write('<h3>C%05d <a href="%s">[KEGG]</a></h3>\n' % (cid, self.kegg().cid2link(cid)))
            self.html_writer.write('Name: %s<br>\n' % self.kegg().cid2name(cid))
            self.html_writer.write('Formula: %s<br>\n' % self.kegg().cid2formula(cid))
            
            pmap = self.cid2pmap_dict[cid]
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0 = pmap.Transform(pH[i], pMg, I[j], T, most_abundant)
                    self.db.Insert('dG0_f', [cid, pH[i], pMg, I[j], T, dG0])
                    self.html_writer.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>f</sub> = %.2f kJ/mol or %.2f kcal/mol<br>\n' % (pH[i], I[j], T, dG0, dG0 / 4.2))
            self.html_writer.write('</p>\n')
        self.db.Commit()
    
    def analyze_all_kegg_reactions(self, pH=[default_pH], pMg=default_pMg, I=[default_I], T=default_T, most_abundant=False):
        self.db.CreateTable('dG0_r', 'rid INT, pH REAL, I REAL, T REAL, dG0 REAL')
        self.db.CreateIndex('dG0_r_idx', 'dG0_r', 'rid, pH, I, T', unique=True)

        self.html_writer.write('<h2><a name=kegg_compounds>&#x394;G<sub>r</sub> of KEGG reactions:</a></h2>')
        for rid in self.kegg().get_all_rids():
            self.html_writer.write('<p>\n')
            self.html_writer.write('<h3>R%05d <a href="%s">[KEGG]</a></h3>\n' % (rid, self.kegg().rid2link(rid)))
            try:
                self.html_writer.write('Definition: %s<br>\n' % self.kegg().rid2reaction(rid).definition)
                self.html_writer.write('Equation:   %s<br>\n' % self.kegg().rid2reaction(rid).equation)
                dG0 = self.estimate_dG_keggrid(rid, pH=pH, I=I, T=T, media=None, most_abundant=most_abundant)
                for i in range(len(pH)):
                    for j in range(len(I)):
                        self.db.Insert('dG0_r', [rid, pH[i], pMg, I[j], T, dG0[i, j]])
                        self.html_writer.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>r</sub> = %.2f kJ/mol<br>\n' % (pH[i], I[j], T, dG0[i, j]))
            except GroupDecompositionError as e:
                self.html_writer.write('Warning, cannot decompose one of the compounds: ' + str(e) + '<br>\n')
            except GroupMissingTrainDataError as e:
                self.html_writer.write('Warning, cannot estimate: ' + str(e) + '<br>\n')
            except KeggParseException as e:
                self.html_writer.write('Warning, cannot parse INCHI: ' + str(e) + '<br>\n')
            except KeyError as e:
                self.html_writer.write('Warning, cannot locate this reaction in KEGG: ' + str(e) + '<br>\n')
            self.html_writer.write('</p>\n')
        self.db.Commit()

    def write_cid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(['cid'] + self.all_group_names)
        for cid in self.kegg().get_all_cids_with_inchi():
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
        for rid in self.kegg().get_all_rids():
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
        self.db.CreateTable('observation', 'cid INT, name TEXT, protons INT, charge INT, nMg INT, dG0_f REAL, use_for TEXT')
        for cid in self.cid2pmap_obs.keys():
            if (cid in self.cid_test_set):
                use_for = 'test'
            else:
                use_for = 'train'
            for (nH, z, nMg, dG0) in self.cid2pmap_obs[cid].ToMatrix():
                self.db.Insert('observation', [cid, self.kegg().cid2name(cid), nH, z, nMg, dG0, use_for])
        
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
        self.cid2pmap_obs = {}
        self.cid_test_set = set()

        for row in self.db.Execute("SELECT * FROM observation"):
            (cid, unused_name, nH, z, mgs, dG0_f, use_for) = row
            self.cid2pmap_obs.setdefault(cid, PseudoisomerMap())
            self.cid2pmap_obs[cid].Add(nH, z, mgs, dG0_f)
            if (use_for == 'test'):
                self.cid_test_set.add(cid)
        
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
        return self.analyze_decomposition(self.kegg().cid2mol(cid))

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
    kegg = Kegg(db)

    if len(sys.argv) < 2:
        html_writer = HtmlWriter('../res/groups.html')
        G = GroupContribution(db=db, html_writer=html_writer, kegg=kegg)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.train("../data/thermodynamics/dG0.csv", use_dG0_format=True)
        G.write_regression_report()
        G.analyze_training_set()
        G.save_cid2pmap()
    else:
        mols = {}
        cid = int(sys.argv[1])
        mols[kegg.cid2name(cid)] = kegg.cid2mol(cid)
        
        G = GroupContribution(db=db, kegg=kegg)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.init()
        #mols['ATP'] = pybel.readstring('smiles', 'C(C1C(C(C(n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O')
        #mols['Tryptophan'] = pybel.readstring('smiles', "c1ccc2c(c1)c(CC(C(=O)O)[NH3+])c[nH]2")
        #mols['Adenine'] = pybel.readstring('smiles', 'c1nc2c([NH2])[n]c[n-]c2n1')
        #mols['Glutamate'] = pybel.readstring('smiles', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OP(=O)([O-])[O-])O)O)O)O')
        #mols['Acetyl-CoA [nH=36]'] = pybel.readstring('smiles', 'CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])(O)=O)n2cnc3c(N)[nH+]cnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')
        #mols['Acetyl-CoA [nH=35]'] = pybel.readstring('smiles', 'CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])(O)=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')
        #mols['Acetyl-CoA [nH=34]'] = pybel.readstring('smiles', 'CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C')

        #mols['acetamide [nH=5]'] = pybel.readstring('smiles', 'CC(=O)N')
        #mols['acetamide [nH=4]'] = pybel.readstring('smiles', 'CC(=O)[NH-]')
        
        #mols['sparteine [nH=27]'] = pybel.readstring('smiles', '[H][C@@]12CCCC[NH+]1C[C@@H]1C[C@H]2C[NH+]2CCCC[C@]12[H]')
        #mols['sparteine [nH=28]'] = pybel.readstring('smiles', '[H][C@@]12CCCC[NH+]1C[C@@H]1C[C@H]2CN2CCCC[C@]12[H]')
        #mols['sparteine [nH=26]'] = pybel.readstring('smiles', '[H][C@@]12CCCCN1C[C@@H]1C[C@H]2CN2CCCC[C@]12[H]')
        #mols['acetyl-CoA a'] = kegg.cid2mol(24)
        #mols['acetyl-CoA b'] = pybel.readstring('smiles', "CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n2cnc3c(N)ncnc23)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C")
        #mols['acetyl-CoA c'] = pybel.readstring('smiles', "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
        #mols['glycylglycine'] = pybel.readstring('smiles', 'C(C(=O)NCC(=O)[O-])[NH3+]')
        #mols['N-Acetylornithine'] = kegg.cid2mol(437)
        #mols['Sinapoyl-CoA'] = kegg.cid2mol(411)
        
        #smarts = pybel.Smarts('C(=O)[N;H1;0]C')
        smarts = None
        
        for key, mol in mols.iteritems():
            #mol.title = key
            print '-'*100
            print key
            if smarts:
                print smarts.findall(mol)
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
            except GroupDecompositionError:
                print "Cannot decompose compound to groups"
            except GroupMissingTrainDataError as e:
                print e.Explain(G)
            mol.draw()
                
