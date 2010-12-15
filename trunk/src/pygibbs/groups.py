#!/usr/bin/python

import sys
import pybel
import openbabel
import csv
import logging
from pylab import nan, find, sum, absolute, array, inv, pinv, plot, log, dot, xlabel, ylabel, figure, hold, text, savefig, zeros, rcParams, mean, std, arange, ylim, xlim, axhspan, axvspan, legend, matrix, exp, contour, meshgrid, clabel
from matplotlib import font_manager
import re
import types

from copy import deepcopy
from toolbox.util import matrixrank, _mkdir
from pygibbs.thermodynamic_constants import R, default_pH, default_pMg, default_I, default_T, default_c0, dG0_f_Mg
from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.feasibility import find_mcmf, LinProgNoSolutionException, find_pCr, thermodynamic_pathway_analysis, pC_to_range
from pygibbs import groups_data, thermodynamic_constants
from pygibbs import group_decomposition
from pygibbs import pseudoisomer
from pygibbs import pseudoisomers_data
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.kegg import KeggParseException, Kegg, parse_bool_field,\
    parse_float_field, parse_vfloat_field, parse_string_field, parse_kegg_file
from toolbox import database, util
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.groups_data import Group
import os

class GroupContributionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        if (type(self.value) == types.StringType):
            return self.value
        else:
            return repr(self.value)
    
class GroupMissingTrainDataError(GroupContributionError):
    def __init__(self, value, missing_groups=[]):
        self.value = value
        self.missing_groups = missing_groups
    def __str__(self):
        if (type(self.value) == types.StringType):
            return self.value
        else:
            return repr(self.value)

class GroupContribution(Thermodynamics):    
    def __init__(self, db, html_writer=None, kegg=None, log_file=sys.stderr):
        Thermodynamics.__init__(self)
        self.db = db

        if html_writer:
            self.FIG_DIR = os.path.basename(html_writer.filename).rsplit('.', 1)[0]
            _mkdir(self.FIG_DIR)
            self.HTML = html_writer
        else:
            self.FIG_DIR = ''
            self.HTML = NullHtmlWriter()

        if kegg:
            self._kegg = kegg
        else:
            self._kegg = Kegg()
        self.bounds = deepcopy(self._kegg.cid2bounds)
        self.hatzi = Hatzi()
        self.hatzi.bounds = self.bounds
            
        self.override_gc_with_measurements = True
        self.cid2pmap_dict = None
            
    def __del__(self):
        self.HTML.close()
    
    def write_gc_tables(self):
        table_names = ["groups", "contribution", "observation"]
        self.HTML.write('<ul>\n')
        for table_name in table_names:
            self.HTML.write('  <li><a href="#%s">Table %s from the database</a></li>\n' % (table_name, table_name))
        self.HTML.write('</ul>\n')
        for table_name in table_names:
            self.HTML.write('<p><h2><a name="%s">Table %s:</a></h2>\n' % (table_name, table_name))
            self.db.Table2HTML(self.HTML, table_name)
            self.HTML.write('</p>\n')
    
    def init(self):
        self.load_groups()
        self.load_contributions()
        self.load_concentrations()
        self.load_training_data()
        self.load_cid2pmap()
        self.kegg()

    def save_cid2pmap(self):
        self.cid2pmap_dict = {}
        logging.info('Calculating the table of chemical formation energies for all KEGG compounds.')
        self.db.CreateTable('gc_cid2prm', 'cid INT, nH INT, z INT, nMg INT, dG0 REAL, estimated BOOL')
        self.db.CreateTable('gc_cid2error', 'id INT, error TEXT')

        for cid in self.kegg().get_all_cids():
            logging.info('Saving KEGG Compound C%05d' % cid)
            
            # If the compound is measured:
            if (cid in self.cid2pmap_obs):
                pmap = self.cid2pmap_obs[cid]
                for (nH, z, nMg, dG0) in pmap.ToMatrix():
                    self.db.Insert('gc_cid2prm', [cid, nH, z, nMg, dG0, False])

            # Try to also estimate the dG0_f using Group Contribution:
            comp = self.kegg().cid2compound(cid)
            if (comp.inchi == None):
                self.db.Insert('gc_cid2error', [cid, 'no InChI exists'])
                continue
            try:
                mol = comp.get_mol() 
            except KeggParseException:
                self.db.Insert('gc_cid2error', [cid, 'cannot determine molecular structure'])
                continue
            try:
                pmap = self.estimate_pmap(mol, ignore_protonations=True)
            except group_decomposition.GroupDecompositionError:
                self.db.Insert('gc_cid2error', [cid, 'cannot decompose into groups'])
                continue
            except GroupMissingTrainDataError:
                self.db.Insert('gc_cid2error', [cid, 'contains groups lacking training data'])
                continue
            self.db.Insert('gc_cid2error', [cid, 'OK'])
            self.cid2pmap_dict[cid] = pmap
            for (nH, z, nMg, dG0) in pmap.ToMatrix():
                self.db.Insert('gc_cid2prm', [cid, int(nH), int(z), int(nMg), dG0, True])
        
        self.db.Commit()
        logging.info('Done writing KEGG compounds.')

    def load_cid2pmap(self):
        self.cid2pmap_dict = {}

        # Now load the data into the cid2pmap_dict:
        for row in self.db.Execute("SELECT cid, nH, z, nMg, dG0 from gc_cid2prm WHERE estimated == 1;"):
            cid, nH, z, nMg, dG0 = row
            self.cid2pmap_dict.setdefault(cid, pseudoisomer.PseudoisomerMap())
            self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)

        cid2pmap_obs = {} # observed formation energies
        for row in self.db.Execute("SELECT cid, nH, z, nMg, dG0 from gc_cid2prm WHERE estimated == 0;"):
            cid, nH, z, nMg, dG0 = row
            cid2pmap_obs.setdefault(cid, pseudoisomer.PseudoisomerMap())
            cid2pmap_obs[cid].Add(nH, z, nMg, dG0)

        # add the observed data to the cid2pmap_dict (and override the estimated data if required)
        for cid in cid2pmap_obs.keys():
            if (self.override_gc_with_measurements or cid not in self.cid2pmap_dict):
                self.cid2pmap_dict[cid] = cid2pmap_obs[cid]
        
    def train(self, pka_fname, obs_fname, use_dG0_format=False):
        l_groupvec_pKa, l_deltaG_pKa, l_names_pKa = self.read_training_data_pKa(pka_fname)
        #l_groupvec_pKa, l_deltaG_pKa, l_names_pKa = [], [], []

        if use_dG0_format:
            l_groupvec_formation, l_deltaG_formation, l_names_formation = self.read_training_data_dG0(obs_fname)
        else:
            l_groupvec_formation, l_deltaG_formation, l_names_formation = self.read_training_data(obs_fname)
        
        l_groupvec = l_groupvec_pKa + l_groupvec_formation
        l_deltaG = l_deltaG_pKa + l_deltaG_formation
        l_names = l_names_pKa + l_names_formation
        self.save_training_data(l_groupvec, l_deltaG, l_names)
        self.group_contributions = self.linear_regression_train()

        logging.info("storing the group contribution data in the database")
        self.db.CreateTable('contribution', 'gid INT, name TEXT, protons INT, charge INT, nMg INT, dG0_gr REAL')
        for i, gc in enumerate(self.group_contributions):
            j = int(self.nonzero_groups[i])
            group = self.groups_data.all_groups[j]
            self.db.Insert('contribution', [j, group.name, group.hydrogens,
                                            group.charge, group.nMg, gc])
            
        self.db.CreateTable('observation', 'cid INT, name TEXT, protons INT, charge INT, nMg INT, dG0_f REAL, use_for TEXT')
        for cid in self.cid2pmap_obs.keys():
            if (cid in self.cid_test_set):
                use_for = 'test'
            else:
                use_for = 'train'
            for (nH, z, nMg, dG0) in self.cid2pmap_obs[cid].ToMatrix():
                self.db.Insert('observation', [cid, self.kegg().cid2name(cid), nH, z, nMg, dG0, use_for])
        
        self.db.Commit()
            
    def load_groups(self, group_fname=None):
        if group_fname:
            self.groups_data = groups_data.GroupsData.FromGroupsFile(group_fname)
            self.groups_data.ToDatabase(self.db)
        else:
            self.groups_data = groups_data.GroupsData.FromDatabase(self.db)
        self.group_decomposer = group_decomposition.GroupDecomposer(self.groups_data)
        self.dissociation = DissociationConstants(self.db, self.HTML,
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
        except group_decomposition.GroupDecompositionError:
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
        
        self.HTML.write('<h2><a name=compounds>List of compounds for training</a></h2>\n')
        self.HTML.write('Source File = %s<br>\n' % obs_fname)
        
        pdata = pseudoisomers_data.PseudoisomersData.FromFile(obs_fname)
        
        counter = 0
        for ps_isomer in pdata:
            name = str(ps_isomer)
            logging.info('Verifying data for %s', name)
            
            self.HTML.write("<h3>%s, %s</h3>\n" % (ps_isomer.name, ps_isomer.ref))

            if not ps_isomer.dG0:
                self.HTML.write('No data for &#x394;G<sub>f</sub><br>\n')
                continue

            if ps_isomer.Skip():
                self.HTML.write('Compound marked as not to be used<br>\n')
                continue
                
            self.HTML.write('&#x394;G<sub>f</sub> = %.2f<br>\n' % ps_isomer.dG0)

            if ps_isomer.cid:
                self.cid2pmap_obs.setdefault(ps_isomer.cid, pseudoisomer.PseudoisomerMap())
                self.cid2pmap_obs[ps_isomer.cid].Add(ps_isomer.hydrogens, ps_isomer.net_charge,
                                                     ps_isomer.magnesiums, ps_isomer.dG0)

            if ps_isomer.Test():
                self.cid_test_set.add(ps_isomer.cid)
                self.HTML.write('Compound marked to be used only for testing (not training)<br>\n')
                continue
            
            elif ps_isomer.Train():
                self.HTML.write('Compound marked to be used for training<br>\n')
            else:
                raise Exception("Unknown usage flag: %s" % ps_isomer.use_for)

            if not ps_isomer.smiles:
                raise Exception("Cannot use compound '%s' for training if it lacks a SMILES string" % ps_isomer.name)
            try:
                self.HTML.write('SMILES = %s<br>\n' % ps_isomer.smiles)
                mol = ps_isomer.Mol()
                mol.removeh()
            except TypeError:
                raise Exception("Invalid smiles: " + ps_isomer.smiles)

            inchi = self.mol2inchi(mol)
            self.HTML.write('INCHI = %s<br>\n' % inchi)
            self.inchi2val_obs[inchi] = ps_isomer.dG0
            mol.title = name
            
            img_fname = self.FIG_DIR + '/train_%05d.png' % counter
            counter += 1
            try:
                self.HTML.embed_molecule_as_png(mol, img_fname)
            except (TypeError, IndexError, AssertionError):
                logging.warning("PyBel cannot draw the compound %s" % name)
                self.HTML.write('WARNING: cannot draw this compound using PyBel\n')

            self.HTML.write('<br>\n')
    
            try:
                decomposition = self.group_decomposer.Decompose(mol, strict=True)
            except group_decomposition.GroupDecompositionError as e:
                logging.error('Cannot decompose one of the compounds in the training set: ' + mol.title)
                raise e
            
            groupvec = decomposition.AsVector()
            l_names.append(name)
            l_groupvec.append(groupvec)
            l_deltaG.append(ps_isomer.dG0)
            self.HTML.write("Decomposition = %s<br>\n" % decomposition)
            
            gc_hydrogens, gc_charge = decomposition.Hydrogens(), decomposition.NetCharge()
            if ps_isomer.hydrogens != gc_hydrogens:
                self.HTML.write("ERROR: Hydrogen count doesn't match: explicit = %d, formula = %d<br>\n" %
                                (ps_isomer.hydrogens, gc_hydrogens))
            if ps_isomer.net_charge != gc_charge:
                self.HTML.write("ERROR: Charge doesn't match: explicit = %d, formula = %d<br>\n" %
                                (ps_isomer.net_charge, gc_charge))
                
        return (l_groupvec, l_deltaG, l_names)        
        
    def read_training_data(self, obs_fname):
        l_groupvec = []
        l_deltaG = []
        l_names = [] # is the list of Molecules (pybel class) used for the regression
        self.inchi2val_obs = {}
        self.cid2pmap_obs = {}
        
        self.HTML.write('<h2><a name=compounds>List of compounds for training</a></h2>\n')
        self.HTML.write('Source File = %s<br>\n' % obs_fname)
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
                self.HTML.write('<h3><a href="%s">C%05d - %s</a></h3>\n' % (self.kegg().cid2link(cid), cid, self.kegg().cid2name(cid)))
                self.HTML.write('Observed value = %f <br>\n' % obs_val)
                
                inchi = self.mol2inchi(mol)
                self.HTML.write('INCHI = %s<br>\n' % inchi)
                self.inchi2val_obs[inchi] = obs_val
                
                groupvec = self.group_decomposer.Decompose(mol).AsVector()
                l_names.append(mol)
                l_groupvec.append(groupvec)
                l_deltaG.append(obs_val)
                self.HTML.write('Decomposition = %s <br>\n' % self.get_decomposition_str(mol))
            except group_decomposition.GroupDecompositionError:
                self.HTML.write('Could not be decomposed<br>\n')
            except KeyError:
                self.HTML.write('Compound has no INCHI in KEGG<br>\n')

        return (l_groupvec, l_deltaG, l_names)        

    def read_training_data_pKa(self, pka_fname):
        
        def smiles2groupvec(smiles, id):
            mol = pybel.readstring('smiles', str(smiles))
            mol.removeh()
            mol.title = id
            self.HTML.embed_molecule_as_png(mol, 'dissociation_constants/%s.png' % id)
            try:
                return self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
            except group_decomposition.GroupDecompositionError as _e:
                logging.warning('Cannot decompose one of the compounds in the training set: %s, %s' % (id, smiles))
                self.HTML.write('%s Could not be decomposed<br>\n' % smiles)
                #raise e
        
        l_groupvec = []
        l_deltaG = []
        l_names = [] # is the list of Molecules' names

        self.dissociation.LoadValuesToDB(pka_fname)
        self.HTML.write('<h2><a name=compounds>List of pKa for training</a></h2>\n')
        self.HTML.write('Source File = %s<br>\n' % pka_fname)
        for cid, step, T, pKa, smiles_below, smiles_above in self.db.Execute("SELECT * FROM pKa WHERE cid IS NOT NULL AND step IS NOT NULL"):
            logging.info("Reading pKa data for C%05d, %s, step %d" % (cid, self.kegg().cid2name(cid), step))
            self.HTML.write('<h3>C%05d - %s - step %d</h3>\n' % (cid, self.kegg().cid2name(cid), step))
            T = T or default_T
            dG0 = R*T*log(10)*pKa
            self.HTML.write('pKa = %.2f, T = %.2f \n</br>' % (pKa, T))
            self.HTML.write('&#x394;G<sub>p</sub> = %.2f<br>\n' % dG0)
            if smiles_below and smiles_above:
                self.HTML.write('SMILES = %s >> %s<br>\n' % (smiles_below, smiles_above))
                decomposition_below = smiles2groupvec(smiles_below, "C%05d_%d_b" % (cid, step))
                decomposition_above = smiles2groupvec(smiles_above, "C%05d_%d_a" % (cid, step))
                if not decomposition_below or not decomposition_above:
                    continue
                groupvec = decomposition_above.AsVector() - decomposition_below.AsVector()
                obs_name = self.kegg().cid2name(cid) + " [%d -> %d]" % \
                    (decomposition_below.NetCharge(), decomposition_above.NetCharge())
                self.HTML.write('<br>\nDecomposition = %s<br>\n' % str(groupvec))
                l_groupvec.append(groupvec)
                l_deltaG.append(dG0)
                l_names.append(obs_name)
            else:
                self.HTML.write('No SMILES provided<br>\n')
        
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
        self.group_matrix = array(X)

        y = []
        for row in self.db.Execute("SELECT * FROM train_observations"):
            y.append(row[0])
        self.obs = array(y)
        
        self.mol_names = []
        for row in self.db.Execute("SELECT name FROM train_molecules"):
            self.mol_names.append(row[0])

    def linear_regression_train(self):
        self.load_training_data()
        
        #### verify that there are no compounds with more than one unique group:
        
        # find all columns that are nonzero only in one row
        unique_groups = find(sum(self.group_matrix != 0, 0) == 1)
        # find rows that have more than one group from that set
        dangerous_rows = find(sum(self.group_matrix[:,unique_groups] != 0, 1) > 1)
        
        # this is dangerous since the weight can be arbitrarily distributed between the two groups
        # and will cause severe miscalculations in predictions of compounds that contain only on of
        # these groups.
        if len(dangerous_rows):
            for j in dangerous_rows:
                logging.warning('The observation %s has more than one unique group, not using it for training' % self.mol_names[j])
            #raise Exception('Cannot perform the regression when there are rows with more than one unique group')
        
        # remove the dangerous rows from the matrix and the obs vector
        self.safe_rows = find(sum(self.group_matrix[:,unique_groups] != 0, 1) < 2)
        self.truncated_group_mat = self.group_matrix[self.safe_rows, :]
        self.truncated_obs = self.obs[self.safe_rows]
        self.truncated_mol_names = [self.mol_names[i] for i in self.safe_rows]
        
        # remove the groups that have no examples from the matrix
        self.nonzero_groups = find(sum(absolute(self.truncated_group_mat), 0) > 0)
        self.truncated_group_mat = self.truncated_group_mat[:, self.nonzero_groups]
        
        inv_corr_mat = pinv(dot(self.truncated_group_mat.T, self.truncated_group_mat))
        #inv_corr_mat = inv(dot(truncated_group_mat.T, truncated_group_mat))
        group_contributions = dot(dot(inv_corr_mat, self.truncated_group_mat.T), self.truncated_obs)
        return group_contributions
    
    def groupvec2val(self, groupvec):
        if (self.group_contributions == None):
            raise Exception("You need to first Train the system before using it to estimate values")

        missing_groups = set(find(groupvec)).difference(set(self.nonzero_groups))
        if (len(missing_groups) == 0):
            groupvec = [groupvec[i] for i in self.nonzero_groups]
            return dot(groupvec, self.group_contributions)
        else:
            raise GroupMissingTrainDataError("can't estimate because some groups have no training data", missing_groups)
    
    def estimate_val(self, mol):
        try:
            groupvec = self.get_groupvec(mol)
        except group_decomposition.GroupDecompositionError as e:
            inchi = self.mol2inchi(mol)
            if (inchi in self.inchi2val_obs):
                return self.inchi2val_obs[inchi]
            else:
                raise e
        
        return self.groupvec2val(groupvec)

    def write_regression_report(self, include_raw_data=False, T=default_T, pH=default_pH):
        if include_raw_data:
            self.db.Table2HTML(self.HTML, 'train_group_matrix')
            self.db.Table2HTML(self.HTML, 'train_observations')
            self.db.Table2HTML(self.HTML, 'train_molecules')
            
        group_matrix_reduced = self.group_matrix[:, self.nonzero_groups]
        self.HTML.write('<h2><a name="regression">Regression</a></h2>\n')
        self.HTML.write('<ul><li>%d compounds</li><li>%d groups</li><li>%d rank</li></ul>\n' % \
                        (group_matrix_reduced.shape[0], group_matrix_reduced.shape[1], matrixrank(group_matrix_reduced)))
        self.HTML.write('<table border="1">\n<tr><td>Name</td><td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td><td>Group Vector</td></tr>')
        for i in range(group_matrix_reduced.shape[0]):
            self.HTML.write('<tr><td>%s</td><td>%.2f</td><td>%s</td></tr>' % \
                            (self.mol_names[i], self.obs[i], ' '.join(['%d' % x for x in group_matrix_reduced[i, :]])))
        self.HTML.write('</table>')
        
        self.HTML.write('<h2><a name="group_contrib">Group Contributions</a></h2>\n')
        self.HTML.write('<table border="1">')
        self.HTML.write('  <tr><td>#</td><td>Group Name</td><td>nH</td><td>charge</td><td>nMg</td><td>&#x394;<sub>gr</sub>G [kJ/mol]</td><td>&#x394;<sub>gr</sub>G\' [kJ/mol]</td><td>Appears in compounds</td></tr>\n')
        for i, group in enumerate(self.nonzero_groups):
            nonzero_group = self.groups_data.all_groups[group]
            dG0_gr = self.group_contributions[i]
            dG0_gr_tag = dG0_gr + R*T*log(10)*pH*nonzero_group.hydrogens
            compound_list_str = ' | '.join([self.mol_names[k] for k in find(group_matrix_reduced[:, i])])
            self.HTML.write('  <tr><td>%d</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%8.2f</td><td>%8.2f</td><td>%s</td></tr>\n' %
                            (i, nonzero_group.name, nonzero_group.hydrogens,
                             nonzero_group.charge, nonzero_group.nMg,
                             dG0_gr, dG0_gr_tag, compound_list_str))
        self.HTML.write('</table>\n')

        self.HTML.write("<p>\nGroups that had no examples in the training set:<br>\n")
        zero_groups = set(range(len(self.groups_data.all_groups))).difference(self.nonzero_groups)
        self.HTML.write("<ol>\n<li>")
        self.HTML.write("</li>\n<li>".join([str(self.groups_data.all_groups[i])
                                            for i in sorted(zero_groups)]))
        self.HTML.write("</li>\n</ol>\n</p>\n")

    def analyze_training_set(self):
        self.write_regression_report()
        
        n_obs = len(self.safe_rows)
        val_obs = []
        val_est = []
        val_err = []
        deviations = []
        
        unique_groups = find(sum(self.truncated_group_mat != 0, 0) == 1)
        row_without_unique_groups = find(sum(self.truncated_group_mat[:, unique_groups] != 0, 1) == 0)
        
        for i in row_without_unique_groups:
            # leave out the row corresponding with observation 'i'
            subset = array(range(0, i) + range((i + 1), n_obs))
            reduced_mat = self.truncated_group_mat[subset, :]
            
            inv_corr_mat = pinv(dot(reduced_mat.T, reduced_mat))
            group_contributions = dot(dot(inv_corr_mat, reduced_mat.T), self.truncated_obs[subset])
            
            estimation = dot(self.truncated_group_mat[i, :], group_contributions)
            error = self.truncated_obs[i] - estimation
            val_obs.append(self.truncated_obs[i])
            val_est.append(estimation)
            val_err.append(error)
            deviations.append((abs(error), self.truncated_mol_names[i], self.truncated_obs[i], estimation, error))
        
        logging.info("writing the table of estimation errors for each compound")
        self.HTML.write('<h2><a name="error_table">Compound Estimation Error</a></h2>\n')
        self.HTML.write('<b>std(error) = %.2f kJ/mol</b>\n' % std(val_err))
        self.HTML.write('<table border="1">')
        self.HTML.write('  <tr><td>Compound Name</td><td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td><td>Error [kJ/mol]</td><td>Remark</td></tr>\n')
        deviations.sort(reverse=True)
        for (_, mol_name, obs, est, err) in deviations:
            self.HTML.write('  <tr><td>%s</td><td>%8.2f</td><td>%8.2f</td><td>%s</td></tr>\n' % (mol_name, obs, err, ""))
        self.HTML.write('</table>\n')
        
        logging.info("Plotting graphs for observed vs. estimated")
        obs_vs_est_fig = figure()
        plot(val_obs, val_est, '.')
        xlabel('Observed (obs)')
        ylabel('Estimated (est)')
        hold(True)
        for (_, mol_name, obs, est, err) in deviations:
            text(obs, est, mol_name, fontsize=4)
        savefig('%s/obs_vs_est.pdf' % self.FIG_DIR, format='pdf')
        self.HTML.write('<h3><a name="obs_vs_est">Observed vs. Estimated</a></h3>\n')
        self.HTML.embed_matplotlib_figure(obs_vs_est_fig, width=1000, height=800)
        
        obs_vs_err_fig = figure()
        plot(val_obs, val_err, '+')
        xlabel('Observed (obs)')
        ylabel('Estimation error (est - obs)')
        hold(True)
        for (_, mol_name, obs, est, err) in deviations:
            text(obs, err, mol_name, fontsize=4)
    
        savefig('%s/obs_vs_err.pdf' % self.FIG_DIR, format='pdf', orientation='landscape')
        self.HTML.write('<h3><a name="obs_vs_err">Observed vs. Error</a></h3>\n')
        self.HTML.embed_matplotlib_figure(obs_vs_err_fig, width=1000, height=800)

    def kegg(self):
        return self._kegg        

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

    def estimate_pmap(self, mol, ignore_protonations=False):
        try:
            all_groupvecs = self.group_decomposer.Decompose(
                mol, ignore_protonations, strict=True).PseudoisomerVectors()
        except group_decomposition.GroupDecompositionError as e:
            raise group_decomposition.GroupDecompositionError(str(e) + "\n" + mol.title + "\n")

        if not all_groupvecs:
            raise group_decomposition.GroupDecompositionError('Found no pseudoisomers for %s'
                                          % mol.title)

        all_missing_groups = []
        pmap = pseudoisomer.PseudoisomerMap()
        for groupvec in all_groupvecs:
            try:
                dG0 = self.groupvec2val(groupvec)
                pmap.AddGroupVector(groupvec, dG0)
            except GroupMissingTrainDataError as e:
                gp_str = ", ".join([self.groups_data.all_group_names[g]
                                    for g in e.missing_groups])
                s = "Species nH = %d, z = %d mg = %d: %s" % (groupvec.Hydrogens(),
                                                             groupvec.NetCharge(),
                                                             groupvec.Magnesiums(),
                                                             gp_str) 
                all_missing_groups.append(s)
        
        if pmap.Empty():
            raise GroupMissingTrainDataError("All species of %s have missing groups:" % mol.title, all_missing_groups)            

        return pmap

    def cid2pmap(self, cid, use_cache=True):
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
            pmap = pseudoisomer.PseudoisomerMap()
            pmap.Add(0, 0, 0, 0)
            return pmap

        if (cid in self.cid2pmap_obs and (self.override_gc_with_measurements or cid in self.cid_test_set)):
            pmap = pseudoisomer.PseudoisomerMap()
            pmap.AddAll(self.cid2pmap_obs[cid])
            return pmap
        else:
            try:
                mol = self.kegg().cid2mol(cid)
            except KeggParseException as e:
                raise MissingCompoundFormationEnergy("Cannot determine molecular structure: " + str(e), cid)
            
            mol.title = "C%05d" % cid
            try:
                return self.estimate_pmap(mol)
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
        for (dG0, unused_nH, z, unused_nMg) in self.cid2pmap(cid).ToMatrix():
            if (z == charge):
                dG0_p0 = dG0
            elif (z == charge + 1):
                dG0_p1 = dG0
        
        if (dG0_p0 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge))
        if (dG0_p1 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge + 1))
        
        return (dG0_p0 - dG0_p1) / (R * T * log(10))

    def estimate_dG0(self, mol, pH=default_pH, pMg=default_pMg,
                     I=default_I, T=default_T):
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).
            
            T in Kelvin
        """
        pmap = self.estimate_pmap(mol)
        return pmap.TransformMatrix(pH, pMg, I, T)
            
    def estimate_dG0_keggcid(self, cid, pH=default_pH, pMg=default_pMg,
                             I=default_I, T=default_T, most_abundant=False): # T = temperature in K
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).            
        """
        pmap = self.cid2pmap(cid)
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
        
        pmaps = [self.cid2pmap(cid) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()

        if type(pH) != types.ListType and type(I) != types.ListType:
            dG0_vector = [pmap.Transform(pH, pMg, I, T, most_abundant) for pmap in pmaps]
            dG0 = dot(stoichiometry_vector, dG0_vector)
            return dG0
        else:
            if type(pH) != types.ListType:
                pH = [float(pH)]
            if type(I) != types.ListType:
                I = [float(I)]
            
            dG0_matrix = zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_vector = [pmap.Transform(pH[i], pMg, I[j], T, most_abundant) for pmap in pmaps]
                    dG0_matrix[i, j] = dot(stoichiometry_vector, dG0_vector)
            
            return dG0_matrix

    def estimate_dG_reaction(self, sparse_reaction, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T, c0=default_c0, media=None, most_abundant=False):
        """
            standard = False means to use the known concentrations of reactants as well
        """
        
        dG0 = self.estimate_dG0_reaction(sparse_reaction, pH, pMg, I, T, most_abundant)
        concentration_vector = [log(self.get_concentration(cid, c0, media)) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()
        
        return dG0 + R * T * dot(stoichiometry_vector, concentration_vector)        

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
        except group_decomposition.GroupDecompositionError:
            raise group_decomposition.GroupDecompositionError("Unable to decompose %s (C%05d) into groups" % (self.kegg().cid2name(cid), cid))
            
    def rid2groupvec(self, rid):
        sparse_reaction = self.kegg().rid2sparse_reaction(rid) 
        group_matrix = matrix([self.cid2groupvec(cid) for cid in sparse_reaction.keys()])
        stoichiometry_vector = matrix(sparse_reaction.values())
        total_groupvec = dot(stoichiometry_vector, group_matrix)
        return total_groupvec.tolist()[0]

    def analyze_all_kegg_compounds(self, pH=[default_pH], pMg=default_pMg, I=[default_I], T=default_T, most_abundant=False):
        self.db.CreateTable('dG0_f', 'cid INT, pH REAL, pMg REAL, I REAL, T REAL, dG0 REAL')
        self.db.CreateIndex('dG0_f_idx', 'dG0_f', 'cid, pH, I, T', unique=True)

        self.HTML.write('<h2><a name=kegg_compounds>&#x394;G<sub>f</sub> of KEGG compounds:</a></h2>')
        for cid in self.cid2pmap_dict.keys():
            self.HTML.write('<p>\n')
            self.HTML.write('<h3>C%05d <a href="%s">[KEGG]</a></h3>\n' % (cid, self.kegg().cid2link(cid)))
            self.HTML.write('Name: %s<br>\n' % self.kegg().cid2name(cid))
            self.HTML.write('Formula: %s<br>\n' % self.kegg().cid2formula(cid))
            
            pmap = self.cid2pmap_dict[cid]
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0 = pmap.Transform(pH[i], pMg, I[j], T, most_abundant)
                    self.db.Insert('dG0_f', [cid, pH[i], pMg, I[j], T, dG0])
                    self.HTML.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>f</sub> = %.2f kJ/mol or %.2f kcal/mol<br>\n' % (pH[i], I[j], T, dG0, dG0 / 4.2))
            self.HTML.write('</p>\n')
        self.db.Commit()
    
    def analyze_all_kegg_reactions(self, pH=[default_pH], pMg=default_pMg, I=[default_I], T=default_T, most_abundant=False):
        self.db.CreateTable('dG0_r', 'rid INT, pH REAL, I REAL, T REAL, dG0 REAL')
        self.db.CreateIndex('dG0_r_idx', 'dG0_r', 'rid, pH, I, T', unique=True)

        self.HTML.write('<h2><a name=kegg_compounds>&#x394;G<sub>r</sub> of KEGG reactions:</a></h2>')
        for rid in self.kegg().get_all_rids():
            self.HTML.write('<p>\n')
            self.HTML.write('<h3>R%05d <a href="%s">[KEGG]</a></h3>\n' % (rid, self.kegg().rid2link(rid)))
            try:
                self.HTML.write('Definition: %s<br>\n' % self.kegg().rid2reaction(rid).definition)
                self.HTML.write('Equation:   %s<br>\n' % self.kegg().rid2reaction(rid).equation)
                dG0 = self.estimate_dG_keggrid(rid, pH=pH, I=I, T=T, media=None, most_abundant=most_abundant)
                for i in range(len(pH)):
                    for j in range(len(I)):
                        self.db.Insert('dG0_r', [rid, pH[i], pMg, I[j], T, dG0[i, j]])
                        self.HTML.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>r</sub> = %.2f kJ/mol<br>\n' % (pH[i], I[j], T, dG0[i, j]))
            except group_decomposition.GroupDecompositionError as e:
                self.HTML.write('Warning, cannot decompose one of the compounds: ' + str(e) + '<br>\n')
            except GroupMissingTrainDataError as e:
                self.HTML.write('Warning, cannot estimate: ' + str(e) + '<br>\n')
            except KeggParseException as e:
                self.HTML.write('Warning, cannot parse INCHI: ' + str(e) + '<br>\n')
            except KeyError as e:
                self.HTML.write('Warning, cannot locate this reaction in KEGG: ' + str(e) + '<br>\n')
            self.HTML.write('</p>\n')
        self.db.Commit()

    def write_cid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(["CID"] + self.all_group_names)
        for cid in self.kegg().get_all_cids_with_inchi():
            try:
                groupvec = self.cid2groupvec(cid)
                csv_file.writerow([cid] + groupvec)
            except group_decomposition.GroupDecompositionError as e:
                print str(e)
                continue
            except KeggParseException as e:
                print str(e)
                continue
            
    def write_rid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(["RID"] + self.all_group_names)
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
            except group_decomposition.GroupDecompositionError as e:
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
        
    def save_contributions(self, filename):
        logging.info("Saving the group contribution data to: %s" % filename)
        csv_output = csv.writer(open(filename, "w"))
        csv_output.writerow(("row type", "ID", "name", "protons", "charge", "Mgs", "dG"))
        for i in range(len(self.group_contributions)):
            j = self.nonzero_groups[i]
            (name, protons, charge, mgs) = self.all_groups[j]
            csv_output.writerow(("CONTRIBUTION", j, name, protons, charge, mgs, self.group_contributions[i]))
            
        for cid in self.cid2pmap_obs:
            for (nH, z, mgs, dG0) in self.cid2pmap_obs[cid].ToMatrix():
                if (cid in self.cid_test_set):
                    use_for = 'test'
                else:
                    use_for = 'train'
                csv_output.writerow(("OBSERVATION", cid, self.kegg().cid2name(cid), nH, z, mgs, dG0, use_for))
            
    def load_contributions(self):
        logging.info("loading the group contribution data from the database")
        self.nonzero_groups = []
        self.group_contributions = []
        self.cid2pmap_obs = {}
        self.cid_test_set = set()
        
        for row in self.db.Execute("SELECT * FROM contribution"):
            (gid, unused_name, unused_protons, unused_charge, unused_mg, dG0_gr) = row
            self.nonzero_groups.append(gid)
            self.group_contributions.append(dG0_gr)
        
        for row in self.db.Execute("SELECT * FROM observation"):
            (cid, unused_name, nH, z, mgs, dG0_f, use_for) = row
            self.cid2pmap_obs.setdefault(cid, pseudoisomer.PseudoisomerMap())
            self.cid2pmap_obs[cid].Add(nH, z, mgs, dG0_f)
            if (use_for == 'test'):
                self.cid_test_set.add(cid)
                
    def read_compound_abundance(self, filename):
        self.db.CreateTable('compound_abundance', 'cid INT, media TEXT, concentration REAL')
        for row in util.ReadCsvWithTitles(filename):
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

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if (len(tokens) == 0):
            return default_value
        if (len(tokens) > 1):
            raise Exception("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    def get_reactions(self, module_name, field_map):
        """
            read the list of reactions from the command file
        """
        
        if ("MODULE" in field_map):
            mid_str = field_map["MODULE"]
            if (mid_str[0] == 'M'):
                mid = int(mid_str[1:])
            else:
                mid = int(mid_str)
            (S, rids, fluxes, cids) = self.kegg().get_module(mid)
            self.HTML.write('<h3>Module <a href=http://www.genome.jp/dbget-bin/www_bget?M%05d>M%05d</a></h3>\n' % (mid, mid))       
        else:
            (S, rids, fluxes, cids) = self.kegg().parse_explicit_module(field_map)

        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        if ("MAP_CID" in field_map):
            for line in field_map["MAP_CID"].split('\t'):
                (cid_before, cid_after) = [int(cid[1:]) for cid in line.split(None, 1)]
                if (cid_before in cids):
                    cids[cids.index(cid_before)] = cid_after
        
        return (S, rids, fluxes, cids)

    def write_metabolic_graph(self, S, rids, cids, svg_fname_graph):
        """
            draw a graph representation of the pathway
        """        
        self.HTML.write('<a href=%s>Graph visualization</a></br>\n' % svg_fname_graph)
        Gdot = self.kegg().draw_pathway(S, rids, cids)
        Gdot.write(svg_fname_graph, prog='dot', format='svg')

    def analyze_profile(self, key, field_map):
        self.HTML.write('<p>\n')
        self.HTML.write('<ul>\n')
        self.HTML.write('<li>Conditions:</br><ol>\n')
        # read the list of conditions from the command file
        conditions = []
        for condition in field_map["CONDITIONS"].split('\t'):
            (media, pH, I, T, c0) = (None, default_pH, default_I, default_T, default_c0)
            media = re.findall("media=([a-zA-Z_]+)", condition)[0]
            if (media == 'None'):
                media = None
            pH = GroupContribution.get_float_parameter(condition, "pH", default_pH)
            I = GroupContribution.get_float_parameter(condition, "I", default_I)
            T = GroupContribution.get_float_parameter(condition, "T", default_T)
            c0 = GroupContribution.get_float_parameter(condition, "c0", default_c0)
            conditions.append((media, pH, I, T, c0))
            self.HTML.write('<li>Conditions: media = %s, pH = %g, I = %g M, T = %g K, c0 = %g</li>\n' % (media, pH, I, T, c0))
        self.HTML.write('</ol></li>\n')
        
        # read the list of methods for calculating the dG
        methods = []
        if (parse_bool_field(field_map, 'MILO', True)):
            methods.append('MILO')
        if (parse_bool_field(field_map, 'HATZI', False)):
            methods.append('HATZI')
        
        # prepare the legend for the profile graph
        legend = []
        dG_profiles = {}
        params_list = []
        for (media, pH, I, T, c0) in conditions:
            for method in methods:
                plot_key = method + ' dG (media=%s,pH=%g,I=%g,T=%g,c0=%g)' % (str(media), pH, I, T, c0)
                legend.append(plot_key)
                dG_profiles[plot_key] = []
                params_list.append((method, media, pH, I, T, c0, plot_key))

        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.kegg().write_reactions_to_html(self.HTML, S, rids, fluxes, cids, show_cids=False)
        self.HTML.write('</ul>')
        self.write_metabolic_graph(self.HTML, S, rids, cids, '%s/%s_graph.svg' % (self.FIG_DIR, key))
        
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = {}
        dG0_r = {}
        dG_f = {}
        dG_r = {}
        for (method, media, pH, I, T, c0, plot_key) in params_list:
            dG0_f[plot_key] = zeros((Nc, 1))
            dG_f[plot_key] = zeros((Nc, 1))
            for c in range(Nc):
                if (method == "MILO"):
                    dG0_f[plot_key][c] = self.cid_to_dG0(cids[c], pH=pH, I=I, T=T)
                elif (method == "HATZI"):
                    dG0_f[plot_key][c] = self.hatzi.cid_to_dG0(cids[c], pH=pH, I=I, T=T)
                else:
                    raise Exception("Unknown dG evaluation method: " + method)
                # add the effect of the concentration on the dG_f (from dG0_f to dG_f)
                dG_f[plot_key][c] = dG0_f[plot_key][c] + R * T * log(self.get_concentration(cids[c], c0, media))
            dG0_r[plot_key] = dot(S, dG0_f[plot_key])
            dG_r[plot_key] = dot(S, dG_f[plot_key])
        
        # plot the profile graph
        rcParams['text.usetex'] = False
        rcParams['legend.fontsize'] = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = 12
        rcParams['lines.linewidth'] = 2
        rcParams['lines.markersize'] = 2
        rcParams['figure.figsize'] = [8.0, 6.0]
        rcParams['figure.dpi'] = 100
        profile_fig = figure()
        hold(True)
        data = zeros((Nr + 1, len(legend)))
        for i in range(len(legend)):
            for r in range(1, Nr + 1):
                data[r, i] = sum(dG_r[legend[i]][:r, 0])
        plot(data)
        legend(legend, loc="lower left")

        for i in range(len(rids)):
            text(i + 0.5, mean(data[i:(i + 2), 0]), rids[i], fontsize=6, horizontalalignment='center', backgroundcolor='white')
        
        xlabel("Reaction no.")
        ylabel("dG [kJ/mol]")
        self.HTML.embed_matplotlib_figure(profile_fig, width=800, heigh=600)
        self.HTML.write('</p>')
    
    def analyze_slack(self, key, field_map):
        self.HTML.write('<p>\n')
        self.HTML.write('<ul>\n')
        self.HTML.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc) 
        c_mid = parse_float_field(field_map, 'C_MID', 1e-3)
        (pH, I, T) = (default_pH, default_I, default_T)
        concentration_bounds = deepcopy(self.kegg().cid2bounds)
        if ("CONDITIONS" in field_map):
            pH = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            I = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            T = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            self.HTML.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                concentration_bounds[cid] = (conc, conc)
                self.HTML.write(', [C%05d] = %g\n' % (cid, conc))
            self.HTML.write('</li>\n')
        self.HTML.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.kegg().write_reactions_to_html(self.HTML, S, rids, fluxes, cids, show_cids=False)
        self.HTML.write('</ul>\n')
        self.write_metabolic_graph(S, rids, cids, '%s/%s_graph.svg' % (self.FIG_DIR, key))
        
        physiological_pC = parse_float_field(field_map, "PHYSIO", 4)
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = zeros((Nc, 1))
        ind_nan = []
        self.HTML.write('<table border="1">\n')
        self.HTML.write('  <td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f [kJ/mol]", "nH", "z"))
        for c in range(Nc):
            cid = cids[c]
            name = self.kegg().cid2name(cid)
            try:
                dG0_f[c] = self.cid_to_dG0(cid, pH, I, T)
                for (nH, z, dG0) in self.cid2pmatrix(cid):
                    self.HTML.write('<tr><td>%05d</td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td>\n' % (cid, name, dG0, nH, z))
            
            except MissingCompoundFormationEnergy:
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = nan
                ind_nan.append(c)
                self.HTML.write('<tr><td>%05d</td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td>\n' % (cid, name))
        self.HTML.write('</table>\n')
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        pC = arange(0, 20, 0.1)
        B_vec = zeros(len(pC))
        #label_vec = [""] * len(pC)
        #limiting_reactions = set()
        for i in xrange(len(pC)):
            c_range = pC_to_range(pC[i], c_mid=c_mid)
            unused_dG_f, unused_concentrations, B = find_mcmf(S, dG0_f, c_range, bounds=bounds)
            B_vec[i] = B
            #curr_limiting_reactions = set(find(abs(dG_r - B) < 1e-9)).difference(limiting_reactions)
            #label_vec[i] = ", ".join(["%d" % rids[r] for r in curr_limiting_reactions]) # all RIDs of reactions that have dG_r = B
            #limiting_reactions |= curr_limiting_reactions

        try:
            unused_dG_f, unused_concentrations, pCr = find_pCr(S, dG0_f, c_mid, bounds=bounds)
        except LinProgNoSolutionException:
            pCr = None
            
        try:
            c_range = pC_to_range(physiological_pC, c_mid=c_mid)
            unused_dG_f, unused_concentrations, B_physiological = find_mcmf(S, dG0_f, c_range, bounds) 
        except LinProgNoSolutionException:
            B_physiological = None

        # plot the profile graph
        rcParams['text.usetex'] = False
        rcParams['legend.fontsize'] = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = 12
        rcParams['lines.linewidth'] = 2
        rcParams['lines.markersize'] = 5
        rcParams['figure.figsize'] = [8.0, 6.0]
        rcParams['figure.dpi'] = 100
        
        slack_fig = figure()
        plot(pC, B_vec, 'b')
        #for i in xrange(len(pC)):
        #    text(pC[i], B_vec[i], label_vec[i], fontsize=6, horizontalalignment='left', backgroundcolor='white')
        xlabel('pC')
        ylabel('slack [kJ/mol]')
        (ymin, _) = ylim()
        (xmin, _) = xlim()
#        broken_barh([(xmin, pCr), (pCr, xmax)], (ymin, 0), facecolors=('yellow', 'green'), alpha=0.3)
        axhspan(ymin, 0, facecolor='b', alpha=0.15)
        title = 'C_mid = %g' % c_mid
        if (pCr != None and pCr < pC.max()):
            title += ', pCr = %.1f' % pCr
            plot([pCr, pCr], [ymin, 0], 'k--')
            text(pCr, 0, 'pCr = %.1f' % pCr, fontsize=8)
            if (pCr < physiological_pC):
                axvspan(pCr, physiological_pC, facecolor='g', alpha=0.3)
        if (B_physiological != None and physiological_pC < pC.max()):
            title += ', slack = %.1f [kJ/mol]' % B_physiological
            plot([xmin, physiological_pC], [B_physiological, B_physiological], 'k--')
            text(physiological_pC, B_physiological, 'B=%.1f' % B_physiological, fontsize=8)
        
        title(title)
        ylim(ymin=ymin)
        self.HTML.embed_matplotlib_figure(slack_fig, width=800, height=600)

        # write a table of the compounds and their dG0_f
        self.HTML.write('<table border="1">\n')
        self.HTML.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f' [kJ/mol]"))
        for c in range(Nc):
            compound = self.kegg().cid2compound(cids[c])
            cid_str = '<a href="%s">C%05d</a>' % (compound.get_link(), compound.cid)
            self.HTML.write('<tr><td>%s</td><td>%s</td><td>%.1f</td>\n' % (cid_str, compound.name, dG0_f[c, 0]))
        self.HTML.write('</table><br>\n')
        
        # write a table of the reactions and their dG0_r
        self.HTML.write('<table border="1">\n')
        self.HTML.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG RID", "Reaction", "flux"))
        for r in range(Nr):
            rid_str = '<a href="http://www.genome.jp/dbget-bin/www_bget?rn:R%05d">R%05d</a>' % (rids[r], rids[r])
            spr = {}
            for c in find(S[r, :]):
                spr[cids[c]] = S[r, c]
            reaction_str = self.kegg().sparse_to_hypertext(spr)
            self.HTML.write('<tr><td>%s</td><td>%s</td><td>%g</td>\n' % (rid_str, reaction_str, fluxes[r]))
        self.HTML.write('</table><br>\n')
        
        self.HTML.write('</p>\n')

    def analyze_margin(self, key, field_map):
        self.HTML.write('<p>\n')
        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.HTML.write('<li>Conditions:</br><ol>\n')
                    
        # The method for how we are going to calculate the dG0
        method = parse_string_field(field_map, "METHOD", "MILO")

        if (method == "MILO"):
            thermodynamics = self
        elif (method == "HATZI"):
            thermodynamics = self.hatzi
        else:
            raise Exception("Unknown dG evaluation method: " + method)
        
        if ("CONDITIONS" in field_map):
            thermodynamics.pH = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            thermodynamics.I = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            thermodynamics.T = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            self.HTML.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (thermodynamics.pH, thermodynamics.I, thermodynamics.T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                thermodynamics.bounds[cid] = (conc, conc)
                self.HTML.write(', [C%05d] = %g\n' % (cid, conc))
            self.HTML.write('</li>\n')
        self.HTML.write('</ol></li>')
        thermodynamics.c_range = parse_vfloat_field(field_map, "C_RANGE", [1e-6, 1e-2])
        thermodynamics.c_mid = parse_float_field(field_map, 'C_MID', 1e-3)
        
        thermodynamic_pathway_analysis(S, rids, fluxes, cids, thermodynamics, self.kegg(), self.HTML)

    def analyze_contour(self, key, field_map):
        pH_list = parse_vfloat_field(field_map, "PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I_list = parse_vfloat_field(field_map, "I", [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        T = parse_float_field(field_map, "T", default_T)
        most_abundant = parse_bool_field(field_map, "ABUNDANT", False)
        formula = parse_string_field(field_map, "REACTION")
        c0 = parse_float_field(field_map, "C0", 1.0)
        media = parse_string_field(field_map, "MEDIA", "None")
        if (media == "None"):
            media = None
            
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        dG_r = self.estimate_dG_reaction(sparse_reaction, pH_list, I_list, T, c0, media, most_abundant)
        contour_fig = figure()
        
        pH_meshlist, I_meshlist = meshgrid(pH_list, I_list)
        CS = contour(pH_meshlist.T, I_meshlist.T, dG_r)       
        clabel(CS, inline=1, fontsize=10)
        xlabel("pH")
        ylabel("Ionic Strength")
        self.HTML.embed_matplotlib_figure(contour_fig, width=800, height=600)
        self.HTML.write('<br>\n' + self.kegg().sparse_to_hypertext(sparse_reaction) + '<br>\n')
        self.HTML.write('</p>')

    def analyze_protonation(self, key, field_map):
        pH_list = parse_vfloat_field(field_map, "PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I = parse_float_field(field_map, "I", default_I)
        T = parse_float_field(field_map, "T", default_T)
        cid = parse_string_field(field_map, "COMPOUND")
        cid = int(cid[1:])
        pmatrix = self.cid2pmatrix(cid)
        data = zeros((len(self.cid), len(pH_list)))
        for j in range(len(pH_list)):
            pH = pH_list[j]
            dG0_array = matrix([-thermodynamic_constants.transform(dG0, nH, z, pH, I, T) / (R * T) \
                                      for (nH, z, dG0) in self.cid])
            dG0_array = dG0_array - max(dG0_array)
            p_array = exp(dG0_array)
            p_array = p_array / sum(p_array)
            data[:, j] = p_array    
        
        protonation_fig = figure()
        plot(pH_list, data.T)
        prop = font_manager.FontProperties(size=10)
        name = self.kegg().cid2name(cid)
        legend(['%s [%d]' % (name, z) for (nH, z, dG0) in pmatrix], prop=prop)
        xlabel("pH")
        ylabel("Pseudoisomer proportion")
        self.HTML.embed_matplotlib_figure(protonation_fig, width=800, height=600)
        self.HTML.write('<table border="1">\n')
        self.HTML.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for (nH, z, dG0) in pmatrix:
            self.HTML.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        self.HTML.write('</table>')
        self.HTML.write('</p>')

    def analyze_pathway(self, filename):
        self.HTML.write("<h1>Pathway analysis using Group Contribution Method</h1>\n")
        self.kegg()
        entry2fields_map = parse_kegg_file(filename)
        
        for key in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[key]
            if (field_map.get("SKIP", "FALSE") == "TRUE"):
                logging.info("skipping pathway: " + key)
                continue
            try:
                self.HTML.write("<b>%s - %s</b>" % (field_map["NAME"], field_map["TYPE"]))
                self.HTML.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (key))
                self.HTML.write('<div id="%s" style="display:none">' % key)
                self.HTML.write('<h2>%s - %s</h2>\n' % (field_map["NAME"], field_map["TYPE"]))
            except KeyError:
                raise Exception("Both the 'NAME' and 'TYPE' fields must be defined for each pathway")

            logging.info("analyzing pathway: " + key)

            if (field_map["TYPE"] == "PROFILE"):     
                self.analyze_profile(key, field_map)
            elif (field_map["TYPE"] == "SLACK"):     
                self.analyze_slack(key, field_map)
            elif (field_map["TYPE"] == "MARGIN"):     
                self.analyze_margin(key, field_map)               
            elif (field_map["TYPE"] == "CONTOUR"):
                self.analyze_contour(key, field_map)
            elif (field_map["TYPE"] == "PROTONATION"):
                self.analyze_protonation(key, field_map)
            else:
                raise Exception("Unknown analysis type: " + field_map["TYPE"])
            self.HTML.write('</div><br>\n')
            
        self.HTML.write('<h4>Measured concentration table:</h4>\n')
        self.HTML.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % ('__concentrations__'))
        self.HTML.write('<div id="%s" style="display:none">' % '__concentrations__')
        self.HTML.write('<p><h2>Abundance</h2>\n')
        self.db.Query2HTML(self.HTML,
                           "SELECT cid, media, 1000*concentration from compound_abundance ORDER BY cid, media",
                           column_names=["CID", "Media", "concentration [mM]"])
        self.HTML.write('</p>\n')
        self.HTML.write('</div><br>\n')

    def add_vectors(self, stoichiometric_vector1, cid_vector1, stoichiometric_vector2, cid_vector2):
        cid_vector = list(set(cid_vector1 + cid_vector2))
        stoichiometric_vector = zeros(len(cid_vector))
        for i in range(len(cid_vector)):
            try:
                i1 = cid_vector1.index(cid_vector[i])
                stoichiometric_vector[i] += stoichiometric_vector1[i1]
            except ValueError:
                pass
            try:
                i2 = cid_vector2.index(cid_vector[i])
                stoichiometric_vector[i] += stoichiometric_vector2[i2]
            except ValueError:
                pass
        
        nonzero_values = find(stoichiometric_vector != 0)
        cid_vector = [cid_vector[i] for i in nonzero_values]
        stoichiometric_vector = stoichiometric_vector[nonzero_values]
        
        return (stoichiometric_vector, cid_vector)
    
    def analyze_decomposition_cid(self, cid):
        return self.analyze_decomposition(self.kegg().cid2mol(cid))

    def analyze_decomposition(self, mol, ignore_protonations=False, strict=False):
        return self.group_decomposer.Decompose(mol, ignore_protonations, strict).ToTableString()
    
    def get_group_contribution(self, name, nH, z, nMg=0):
        gr = Group(None, name, nH, z, nMg)
        index_full = self.groups_data.Index(gr)
        index_short = self.nonzero_groups.index(index_full)
        dG0_gr = self.group_contributions[index_short]
        return dG0_gr
    
    def calc_pKa_group(self, name, nH, z, nMg=0):
        dG0_gr0 = self.get_group_contribution(name, nH-1, z-1, nMg)
        dG0_gr1 = self.get_group_contribution(name, nH, z, nMg)
        pKa = (dG0_gr0 - dG0_gr1)/(R*default_T*log(10))
        return pKa
    
    def calc_pK_Mg_group(self, name, nH, z, nMg):
        dG0_gr0 = self.get_group_contribution(name, nH, z-2, nMg-1)
        dG0_gr1 = self.get_group_contribution(name, nH, z, nMg)
        pK_Mg = (dG0_gr0 + dG0_f_Mg - dG0_gr1)/(R*default_T*log(10))
        return pK_Mg
    
#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
    
if __name__ == '__main__':
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    if True:
        html_writer = HtmlWriter('../res/dG0_train.html')
        G = GroupContribution(db, html_writer)
        G.load_groups("../data/thermodynamics/groups_species.csv")
        G.train('../data/thermodynamics/pKa_with_cids.csv', "../data/thermodynamics/dG0.csv", use_dG0_format=True)
        G.analyze_training_set()
        G.save_cid2pmap()
    else:
        html_writer = HtmlWriter('../res/dG0_test.html')
        G = GroupContribution(db, html_writer)
        G.init()
        G.load_groups("../data/thermodynamics/groups_species.csv")
        
        mols = {}
        #mols['ATP'] = pybel.readstring('smiles', 'C(C1C(C(C(n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O')
        #mols['Tryptophan'] = pybel.readstring('smiles', "c1ccc2c(c1)c(CC(C(=O)O)[NH3+])c[nH]2")
        #mols['Adenine'] = pybel.readstring('smiles', 'c1nc2c([NH2])[n]c[n-]c2n1')
        #mols['Glutamate'] = pybel.readstring('smiles', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OP(=O)([O-])[O-])O)O)O)O')
        mols['Tryptophan [0]'] = pybel.readstring('smiles', 'c1ccc2c(c1)c(CC(C(=O)[O-])[NH3+])c[nH]2')
        mols['Tryptophan [-1]'] = pybel.readstring('smiles', 'c1ccc2c(c1)c(CC(C(=O)[O-])N)c[nH]2')
        
        smarts = pybel.Smarts("[n;H1;X3;+1]")
        
        for key, mol in mols.iteritems():
            mol.title = key
            #mol.draw()
            print '-'*100
            print key
            print smarts.findall(mol)
            print G.analyze_decomposition(mol)
            #pmap = G.estimate_pmap(mol, ignore_protonations=False)
            #print pmap

