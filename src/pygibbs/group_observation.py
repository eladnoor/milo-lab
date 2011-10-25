import logging
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.group_vector import GroupVector
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T,\
    default_pMg
import csv
from pygibbs.nist import Nist
import json
from pygibbs.dissociation_constants import MissingDissociationConstantError
from toolbox.linear_regression import LinearRegression

class GroupObservation(object):
    
    def __init__(self):
        self.sparse = None
        self.dG0 = None
        self.name = None
        self.kegg_id = None
        self.obs_type = None # can be 'formation', 'reaction', 'acid-base' or 'Mg'

    def NetCharge(self):
        return self.groupvec.NetCharge()

    def Hydrogens(self):
        return self.groupvec.Hydrogens()

    def Magnesiums(self):
        return self.groupvec.Magnesiums()

    def ToDatabase(self, db, table_name):
        db.Insert(table_name, self.ToCsvRow())
    
    @staticmethod 
    def FromDatabaseRow(row):
        obs = GroupObservation()
        obs.name = row['name']
        obs.kegg_id = row['kegg_id']
        obs.obs_type = row['obs_type']
        obs.dG0 = row['dG0']
        obs.sparse = json.loads(row['stoichiometry'])
        return obs
        
    def ToCsvRow(self):
        return [self.name, self.kegg_id, self.obs_type, self.dG0, 
                json.dumps(self.sparse)]
    
    def __hash__(self):
        return hash(json.dumps(self.sparse))
    
    def __eq__(self, other):
        return self.sparse == other.sparse
    
    def Normalize(self):
        """
            In self.sparse, the first non-zero coefficient must be 1,
            where the order is defined alphabetically by the keys.
        """
        if self.sparse == {}:
            return
        factor = 1.0 / self.sparse[min(self.sparse.keys())]
        self.sparse = dict([(k,v*factor) for (k,v) in self.sparse.iteritems()])
        self.dG0 *= factor
        
class GroupObervationCollection(object):
    
    def __init__(self, db, html_writer, group_decomposer, dissociation,
                 transformed=False):
        self.kegg = Kegg.getInstance()
        self.db = db
        self.groups_data = group_decomposer.groups_data
        self.group_decomposer = group_decomposer
        self.observations = []
        self.id2gv = {}
        self.n_groups = len(self.groups_data.all_groups)
        self.groupvec_titles = ["g%d" % i for i in xrange(self.n_groups)]
        self.dissociation = dissociation
        self.transformed = transformed
        self.use_pkas = False
        self.html_writer = html_writer
        if self.transformed:
            self.TRAIN_GROUPS_TABLE_NAME = 'bgc_train_groups'
            self.PSEUDOISOMER_GV_TABLE_NAME = 'bgc_pseudoisomer_groupvector'
            self.OBSERVATION_TABLE_NAME = 'bgc_observations'
        else:
            self.TRAIN_GROUPS_TABLE_NAME = 'pgc_train_groups'
            self.PSEUDOISOMER_GV_TABLE_NAME = 'pgc_pseudoisomer_groupvector'
            self.OBSERVATION_TABLE_NAME = 'pgc_observations'
            
    def FromFiles(self):
        if self.use_pkas:
            self.html_writer.write('<br><b>List of pKa for training</b>')
            self.html_writer.insert_toggle(start_here=True)
            self.AddDissociationTable()
            self.html_writer.div_end()

        self.html_writer.write('</br><b>List of NIST reactions for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.AddNistDatabase()
        self.html_writer.div_end()
        
        self.html_writer.write('<br><b>List of compounds for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.AddFormationEnergies()
        self.html_writer.div_end()
    
    def AddPseudoisomer(self, pseudoisomer_id, mol):
        """
            Decomposes a molecule and calculates its GroupVector.
            Adds the GroupVector to the dictionary according to the given compound_id.
            
            Returns:
                True only if successful
        """
        if mol is None:
            return False
        
        try:
            mol.RemoveHydrogens()
            decomposition = self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
        except GroupDecompositionError as _e:
            logging.warning('Cannot decompose one of the compounds in the '
                            'training set: %s, %s' % (id, str(mol)))
            self.html_writer.write('%s Could not be decomposed</br>\n' % str(mol))
            return False
        groupvec = decomposition.AsVector()
        self.id2gv[pseudoisomer_id] = groupvec
        return True
    
    def AddTrainingExample(self, sparse, dG0, name, kegg_id, obs_type):
        obs = GroupObservation()
        obs.sparse = sparse
        obs.dG0 = dG0
        obs.name = name
        obs.kegg_id = kegg_id
        obs.obs_type = obs_type
        obs.Normalize()
        self.observations.append(obs)

    def AddDissociationTable(self):
        if self.transformed:
            raise Exception("Cannot add pKas to the training set when running"
                            " in 'transformed' mode")
        for cid in self.dissociation.GetAllCids():
            diss_table = self.dissociation.GetDissociationTable(cid)
            if diss_table is None:
                continue
            name = self.kegg.cid2name(cid)
            logging.debug("Reading pKa data for C%05d, %s:" % (cid, name))
            self.html_writer.write('<b>C%05d - %s</b>\n' % (cid, name))
            for key in diss_table:
                nH_above, nH_below, nMg_above, nMg_below = key
                id_below = "C%05d_H%d_Mg%d" % (cid, nH_below, nMg_below)
                id_above = "C%05d_H%d_Mg%d" % (cid, nH_above, nMg_above)
                sparse = {id_below:1, id_above:-1}
                
                ddG0, _ref = diss_table.ddGs[key]
                mol_above = diss_table.GetMol(nH_above, nMg_above)
                mol_below = diss_table.GetMol(nH_below, nMg_below)
                
                if not self.AddPseudoisomer(id_below, mol_below):
                    continue
                if not self.AddPseudoisomer(id_above, mol_above):
                    continue
        
                if nH_above != nH_below:
                    s = "nH = %d -> %d" % (nH_below, nH_above)
                elif nMg_above != nMg_below:
                    s = "nMg = %d -> %d" % (nMg_below, nMg_above)
                logging.debug("\t" + s)
                html_text = ""
                html_text += '<font size="1">\n'
                html_text += '&#x394;&#x394;G = %.2f, %s</br>\n' % (ddG0, s)
                html_text += 'SMILES = %s >> %s</br>\n' % (mol_above.ToSmiles(), mol_below.ToSmiles())
                html_text += '%s >> %s</br>\n' % (mol_above.ToSVG(), mol_below.ToSVG())
                html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
                html_text += '</font>\n'
                self.html_writer.write(html_text)

                if nH_above != nH_below:
                    if nH_above != nH_below-1:
                        raise Exception("a pKa must represent a difference of exactly 1 hydrogen")
                    fullname = '%s [nH=%d->%d]' % (name, nH_above, nH_below)
                    self.AddTrainingExample(sparse, ddG0, fullname, kegg_id="C%05d" % cid, obs_type='acid-base')
                elif nMg_above != nMg_below:
                    if nMg_above != nMg_below-1:
                        raise Exception("a pK_Mg must represent a difference of exactly 1 magnesium")
                    fullname = '%s [nMg=%d->%d]' % (name, nMg_above, nMg_below)
                    raise Exception("Mg dissociation constants are not supported in PGC")
                    #self.AddTrainingExample(sparse, ddG0, fullname, cid, obs_type='Mg')
                    # TODO: solve the problems with Mg species !!!
                    # there seems to be a problem in decomposition of ADP:Mg (and maybe other forms that have Mg ions)

    def AddBiochemicalFormationEnergies(self, 
            obs_fname="../data/thermodynamics/formation_energies_transformed.csv"):
        for row in csv.DictReader(open(obs_fname, 'r')):
            if row['use for'] != 'training':
                continue
            cid = int(row['cid'])
            dG0_prime = float(row["dG'0"])
            pH = float(row['pH'])
            I = float(row['I'])
            pMg = float(row['pMg'])
            T = float(row['T'])
            ref = row['compound_ref']
            kegg_id = "C%05d" % cid
            
            name = "%s, %s" % (self.kegg.cid2name(cid), ref)
            logging.debug('Adding the formation energy of %s', name)
            mol = self.dissociation.GetAnyMol(cid)
            if not self.AddPseudoisomer(kegg_id, mol):
                raise Exception('Cannot decompose %s from the training set: %s'
                              % (kegg_id, mol.ToSmiles()))
    
            sparse = {kegg_id:1}
            self.AddTrainingExample(sparse, dG0_prime, name, kegg_id=kegg_id, obs_type='formation')

            html_text = ""
            html_text += "<b>%s (%s), %s</b></br>\n" % (name, kegg_id, ref)
            html_text += '<font size="1">\n'
            html_text += '&#x394;<sub>f</sub> G\'<sup>0</sup> = %.2f, ' % dG0_prime
            html_text += 'pH = %g, I = %g, pMg = %g, T = %g</br>\n' % (pH, I, pMg, T)
            html_text += 'SMILES = %s</br>\n' % (mol.ToSmiles())
            html_text += '%s</br>\n' % (mol.ToSVG())
            html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
            html_text += '</font>\n'
            self.html_writer.write(html_text)
    
    def AddFormationEnergies(self,
            obs_fname="../data/thermodynamics/formation_energies.csv",
            pH=default_pH, I=default_I, pMg=default_pMg, T=default_T):
        """
            Add observations based on a table of derived chemical formation energies.
            If working in self.transformed mode, the dG0_f is Legendre-transformed
            to dG0'_f using the table of dissociation constants.
        """
        
        train_species = PsuedoisomerTableThermodynamics.FromCsvFile(obs_fname, label='training')
        for cid in train_species.get_all_cids():
            pmap = train_species.cid2PseudoisomerMap(cid)
            pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
            if len(pmatrix) != 1:
                raise Exception("multiple training species for C%05d" % cid)
            nH, charge, nMg, dG0 = pmatrix[0]
            name = "%s (%d)" % (self.kegg.cid2name(cid), nH)
            logging.debug('Adding the formation energy of %s', name)
            diss_table = self.dissociation.GetDissociationTable(cid, 
                                                        create_if_missing=True)
            if diss_table is None:
                raise Exception("%s [C%05d, nH=%d, nMg=%d] does not have a " 
                                "dissociation table"
                                % (name, cid, nH, nMg))

            diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
            dG0_prime = diss_table.Transform(pH, I, pMg, T)
            mol = diss_table.GetMol(nH=nH, nMg=nMg)
            if mol is None:
                raise Exception("%s [C%05d, nH=%d, nMg=%d] does not have a SMILES "
                                "expression in the dissociation constant table" 
                                % (name, cid, nH, nMg))
            mol.SetTitle(name)
            if self.transformed:
                pseudoisomer_id = "C%05d" % (cid)
            else:
                pseudoisomer_id = "C%05d_H%d_Mg%d" % (cid, nH, nMg)
            if not self.AddPseudoisomer(pseudoisomer_id, mol):
                raise Exception('Cannot decompose %s from the training set: %s'
                              % (pseudoisomer_id, mol.ToSmiles()))
            sparse = {pseudoisomer_id:1}
    
            html_text = ""
            html_text += "<b>%s (C%05d), %s</b></br>\n" % \
                        (name, cid, train_species.cid2SourceString(cid))
            html_text += '<font size="1">\n'
            if self.transformed:
                html_text += '&#x394;<sub>f</sub> G\'<sup>0</sup> = %.2f, ' % dG0_prime
                html_text += 'pH = %g, I = %g, pMg = %g, T = %g</br>\n' % (pH, I, pMg, T)
            else:
                html_text += '&#x394;<sub>f</sub> G<sup>0</sup> = %.2f, ' % dG0
                html_text += 'nH = %d, nMg = %d</br>\n' % (nH, nMg)
            html_text += 'SMILES = %s</br>\n' % (mol.ToSmiles())
            html_text += '%s</br>\n' % (mol.ToSVG())
            html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
    
            if self.transformed:
                self.AddTrainingExample(sparse, dG0_prime, name, kegg_id="C%05d" % cid, obs_type='formation')
            else:
                self.AddTrainingExample(sparse, dG0, name, kegg_id="C%05d" % cid, obs_type='formation')
                
                # check that the nH and z of the decomposition matches
                # the indicated values from the CSV file
                gc_nH = self.id2gv[pseudoisomer_id].Hydrogens()
                gc_charge = self.id2gv[pseudoisomer_id].NetCharge()
                if nH != gc_nH:
                    s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, decomposition = %d' % (
                        nH, gc_nH)
                    logging.error(s)
                    html_text += s + '</br>\n'
                if charge != gc_charge:
                    s = 'ERROR: Charge doesn\'t match: explicit = %d, decomposition = %d' % (
                        charge, gc_charge)
                    logging.error(s)
                    html_text += s + '</br>\n'
            html_text += '</font>\n'
            self.html_writer.write(html_text)
    
    def AddNistDatabase(self):
        """
            Add the observations based on equilibrium constants from the NIST database.
            If using non-transformed group contribution, it is required to reverse Legendre-
            transform the data in order to get chemical reaction energies.
            
            This methods tries to use the same pseudoisomers for each compound as
            the table of formation energies, in order to have less columns in the 
            final stoichiometric matrix.
        """
        # collect the formation energies of the anchored compounds (i.e. 'testing')
        # in a dictionary according to CIDs.
        cid2dG0 = {} # this is dG0' if self.transformed, or dG0 if not.
        if self.transformed:
            obs_fname = '../data/thermodynamics/formation_energies_transformed.csv'
        else:
            obs_fname = '../data/thermodynamics/formation_energies.csv'
            cid2nH_nMg = self.dissociation.GetCid2nH_nMg(pH=7, I=0, pMg=14, T=default_T)
        
        for row in csv.DictReader(open(obs_fname, 'r')):
            cid = int(row['cid'])
            if row['use for'] == 'testing':
                if self.transformed:
                    cid2dG0[cid] = float(row["dG'0"])
                else:
                    cid2dG0[cid] = float(row["dG0"])
                    cid2nH_nMg[cid] = (int(row["nH"]), int(row["nMg"])) # override defaults

        if not self.transformed:
            # override the nH and nMg of the pseudoisomers used for training
            # so it will be consistent with the formation energy table
            train_species = PsuedoisomerTableThermodynamics.FromCsvFile(
                '../data/thermodynamics/formation_energies.csv', label='training')
            for cid in train_species.get_all_cids():
                pmap = train_species.cid2PseudoisomerMap(cid)
                pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
                if len(pmatrix) != 1:
                    raise Exception("C%05d has more than one species in the training set" % cid)
                nH, _z, nMg, _dG0 = pmatrix[0]
                cid2nH_nMg[cid] = (nH, nMg) # override defaults
                    
        nist = Nist()
        
        # add all the non-anchored compounds to the dictionary of GroupVectors
        for cid in set(nist.GetAllCids()).difference(cid2dG0.keys()):
            name = self.kegg.cid2name(cid)
            if self.transformed:
                mol = self.dissociation.GetAnyMol(cid)
                pseudoisomer_id = "C%05d" % cid
            else:
                nH, nMg = cid2nH_nMg.get(cid, (0,0))
                mol = self.dissociation.GetMol(cid, nH, nMg)
                pseudoisomer_id = "C%05d_nH%d_nMg%d" % (cid, nH, nMg)
            
            if mol is None:
                logging.warning('C%05d has no explicit formula' % cid)
                continue
                
            if not self.AddPseudoisomer(pseudoisomer_id, mol):
                raise Exception('Cannot decompose %s from the training set: %s'
                              % (pseudoisomer_id, mol.ToSmiles()))
            
            html_text = ""
            html_text += "<b>%s (%s)</b></br>\n" % (name, pseudoisomer_id)
            html_text += '<font size="1">\n'
            html_text += 'SMILES = %s</br>\n' % (mol.ToSmiles())
            html_text += '%s</br>\n' % (mol.ToSVG())
            html_text += 'Stoichiometry = %s</br>\n' % str({pseudoisomer_id:1})
            html_text += '</font>\n'
            self.html_writer.write(html_text)

        # create a dictionary from each unique reaction to the list of measured dG0'
        # and subtract from dG0' the formation energies of the anchored compounds
        for r, nist_row_data in enumerate(nist.SelectRowsFromNist()):
            name="NIST%03d" % r
            dG0 = nist_row_data.dG0_r
            if not self.transformed:
                try:
                    dG0 -= self.dissociation.ReverseTransformNistRow(
                                                    nist_row_data, cid2nH_nMg)
                except MissingDissociationConstantError as e:
                    logging.warning('Cannot reverse transform NIST%03d because of'
                    ' of a missing dissociation constant for C%05d' % (r, e.cid))
                    continue
            
            sparse = {}
            for cid, coeff in nist_row_data.reaction.iteritems():
                if cid in cid2dG0:
                    dG0 -= cid2dG0[cid] * coeff # subtract the effect of the anchored compounds
                else:
                    if self.transformed:
                        pseudoisomer_id = "C%05d" % cid
                    else:
                        nH, nMg = cid2nH_nMg[cid]
                        pseudoisomer_id = "C%05d_nH%d_nMg%d" % (cid, nH, nMg)
                    sparse[pseudoisomer_id] = coeff

            if not set(sparse.keys()).issubset(self.id2gv.keys()):
                logging.warning('Cannot reverse transform NIST%03d because it'
                                ' involves implicit-formula reactants' % (r, e.cid))
                continue
            
            self.AddTrainingExample(sparse, dG0, name=name,
                                    kegg_id="", obs_type="reaction")
            
            html_text = ""
            html_text += "<b>%s</b></br>\n" % name
            html_text += '<font size="1">\n'
            if self.transformed:
                symbol = "&#x394;<sub>r</sub> G'<sup>0</sup>"
            else:
                symbol = "&#x394;<sub>r</sub> G<sup>0</sup>"
            html_text += '%s = %.2f, pH = %g, I = %g, pMg = %g, T = %g</br>\n' % \
                         (symbol, nist_row_data.dG0_r, nist_row_data.pH,
                         nist_row_data.I, nist_row_data.pMg, nist_row_data.T)
            html_text += "Original reaction = %s</br>\n" + nist_row_data.reaction.to_hypertext()
            html_text += '%s (truncated) = %.2f</br>\n' % (symbol, dG0)
            html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
            html_text += '</font>\n'
            self.html_writer.write(html_text)

    def ToDatabase(self):
        # This table is used only as an output for checking results
        # it is not used elsewhere in the code
        self.db.CreateTable(self.TRAIN_GROUPS_TABLE_NAME, 
                            'name TEXT, nH INT, z INT, nMg INT',
                            drop_if_exists=True)
        for group in self.groups_data.all_groups:
            self.db.Insert(self.TRAIN_GROUPS_TABLE_NAME,
                    [group.name, group.hydrogens, group.charge, group.nMg])

        # the PSEUDOISOMER_GV_TABLE_NAME will contain a mapping from the pseudoisomer
        # IDs to the corresponding GroupVectors
        self.db.CreateTable(self.PSEUDOISOMER_GV_TABLE_NAME,
                            'id TEXT, groupvec TEXT',
                            drop_if_exists=True)
        for id, gv in self.id2gv.iteritems():
            self.db.Insert(self.PSEUDOISOMER_GV_TABLE_NAME, [id, gv.ToJSONString()])

        # the table 'group_observation' will contain all the observed data
        # that is used for training later
        titles = ['name TEXT', 'kegg_id TEXT', 'obs_type TEXT', 'dG0 REAL', 'stoichiometry TEXT']
        self.db.CreateTable(self.OBSERVATION_TABLE_NAME, ','.join(titles),
                            drop_if_exists=True)

        for obs in self.observations:
            obs.ToDatabase(self.db, table_name=self.OBSERVATION_TABLE_NAME)
        self.db.Commit()
        
    def ToCSV(self, gv_fname, obs_fname):
        # write the mapping from IDs to group vectors
        csv_writer = csv.writer(open(gv_fname, 'w'))
        csv_writer.writerow(['id', 'groupvec'])
        for id, gv in self.id2gv.iteritems():
            csv_writer.writerow([id, gv.ToJSONString()])

        # write all observations
        csv_writer = csv.writer(open(obs_fname, 'w'))
        csv_writer.writerow(['name', 'kegg_id', 'obs_type', 'dG0', 'stoichiometry'])
        for obs in self.observations:
            csv_writer.writerow(obs.ToCsvRow())
            
    def FromDatabase(self):
        self.id2gv = {}
        for row in self.db.DictReader(self.PSEUDOISOMER_GV_TABLE_NAME):
            gv = GroupVector.FromJSONString(self.groups_data, row['groupvec'])
            self.id2gv[row['id']] = gv
        
        self.observations = []
        for row in self.db.DictReader(self.OBSERVATION_TABLE_NAME):
            obs = GroupObservation.FromDatabaseRow(row)
            self.observations.append(obs)
        
    def GetRegressionData(self):
        """ 
            Returns the regression matrix and corresponding observation data 
            and names as a tuple: (A, b, names)
            Next step is to solve Ax = b
            
            Makes sure that there are no duplicates (two identical rows in A)
            by averaging the observed values.
        """

        # first create the group matrix G (rows=pseudoisomer, cols=groups)
        id_list = sorted(self.id2gv.keys())
        id2index = dict([(id, i) for (i, id) in enumerate(id_list)])
        
        G = np.zeros((len(id_list), len(self.groups_data.all_groups)))
        for i_pseudoisomer, id in enumerate(id_list):
            for i_group, coeff in enumerate(self.id2gv[id]):
                G[i_pseudoisomer, i_group] = coeff
                
        # then create the stoichiometric matrix S (rows=observation, cols=pseudoisomers)
        S = np.zeros((len(self.observations), len(id_list)))
        dG_vec = np.zeros((len(self.observations), 1))
        for i_observation, obs in enumerate(self.observations):
            for id, coeff in obs.sparse.iteritems():
                i_pseudoisomer = id2index[id]
                S[i_observation, i_pseudoisomer] = coeff
            dG_vec[i_observation, 0] = obs.dG0
            
        # "squeeze" S and dG_vec such that repeating rows will be combined into a single
        # row, and its observation will be their mean dG0.
        S, dG_vec = LinearRegression.RowUnique(S, dG_vec)
        
        obs_types = [obs.obs_type for obs in self.observations]
        names = [obs.name for obs in self.observations]

        regression_matrix = np.dot(S, G)
        return regression_matrix, dG_vec, obs_types, names
