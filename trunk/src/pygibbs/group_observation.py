import logging
import os
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.group_vector import GroupVector
from pygibbs.nist_regression import NistRegression
from toolbox.molecule import Molecule
from toolbox.html_writer import NullHtmlWriter
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T

class GroupObservation(object):
    
    def __init__(self):
        self.groupvec = None
        self.dG0 = None
        self.name = None
        self.id = None
        self.obs_type = None # can be 'formation', 'reaction', 'acid-base' or 'Mg'

    def NetCharge(self):
        return self.groupvec.NetCharge()

    def Hydrogens(self):
        return self.groupvec.Hydrogens()

    def Magnesiums(self):
        return self.groupvec.Magnesiums()

    def ToDatabase(self, db):
        db.Insert('group_observations', list(self.groupvec) +
                  [self.dG0, self.name, self.id, self.obs_type])
        
    def __hash__(self):
        return hash(str(self.groupvec))
    
    def __eq__(self, other):
        return self.groupvec == other.groupvec
        
class GroupObervationCollection(object):
    
    def __init__(self, db, html_writer, group_decomposer):
        self.kegg = Kegg.getInstance()
        self.db = db
        self.groups_data = group_decomposer.groups_data
        self.group_decomposer = group_decomposer
        self.observations = []
        self.n_groups = len(self.groups_data.all_groups)
        self.groupvec_titles = ["g%d" % i for i in xrange(self.n_groups)]

        self.html_writer = html_writer
        if html_writer and html_writer.filename:
            self.FIG_DIR = os.path.basename(html_writer.filename).rsplit('.', 1)[0]
        else:
            self.FIG_DIR = ''
        
    def Add(self, groupvec, dG0, name, id, obs_type):
        obs = GroupObservation()
        obs.groupvec = groupvec
        obs.dG0 = dG0
        obs.name = name
        obs.id = id
        obs.obs_type = obs_type
        self.observations.append(obs)

    def SmilesToGroupvec(self, smiles, id):
        mol = Molecule.FromSmiles(str(smiles))
        mol.RemoveHydrogens()
        mol.SetTitle(id)
        
        #try:
        self.html_writer.write(mol.ToSVG())
        #except (TypeError, AssertionError): # The Mg ions cannot be drawn by OASA 
        #    pass
        
        try:
            return self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
        except GroupDecompositionError as _e:
            logging.warning('Cannot decompose one of the compounds in the training set: %s, %s' % (id, smiles))
            self.html_writer.write('%s Could not be decomposed</br>\n' % smiles)

    def AddDissociationTable(self):
        dissociation = DissociationConstants.FromFile()
        for cid in dissociation.GetAllCids():
            diss_table =  dissociation.GetDissociationTable(cid)
            name = self.kegg.cid2name(cid)
            logging.info("Reading pKa data for C%05d, %s:" % (cid, name))
            self.html_writer.write('<h3>C%05d - %s</h3>\n' % (cid, name))
            for key in diss_table:
                nH_above, nH_below, nMg_above, nMg_below = key
                ddG0, _ref = diss_table.ddGs[key]
                smiles_above = diss_table.smiles_dict[nH_above, nMg_above]
                smiles_below = diss_table.smiles_dict[nH_below, nMg_below]
                
                if not smiles_above or not smiles_below:
                    continue
        
                if nH_above != nH_below:
                    s = "nH = %d -> %d" % (nH_below, nH_above)
                elif nMg_above != nMg_below:
                    s = "nMg = %d -> %d" % (nMg_below, nMg_above)
                logging.info("\t" + s)
                self.html_writer.write('%s: &#x394;G = %.2f</br>\n' 
                                       % (s, ddG0))
                self.html_writer.write('SMILES = %s >> %s</br>\n' % (smiles_below, smiles_above))
                decomposition_below = self.SmilesToGroupvec(smiles_below, 
                                "C%05d_b_H%d_Mg%d" % (cid, nH_below, nMg_below))
                decomposition_above = self.SmilesToGroupvec(smiles_above, 
                                "C%05d_a_H%d_Mg%d" % (cid, nH_above, nMg_above))
                if not decomposition_below or not decomposition_above:
                    continue
                groupvec = decomposition_below.AsVector() - decomposition_above.AsVector()
                self.html_writer.write('</br>\nDecomposition = %s</br>\n' % str(groupvec))
                
                if nH_above != nH_below:
                    if nH_above != nH_below-1:
                        raise Exception("a pKa must represent a difference of exactly 1 hydrogen")
                    fullname = '%s [nH=%d]' % (name, nH_above)
                    self.Add(groupvec, ddG0, fullname, "C%05d" % cid, 'acid-base')
                elif nMg_above != nMg_below:
                    if nMg_above != nMg_below-1:
                        raise Exception("a pK_Mg must represent a difference of exactly 1 magnesium")
                    fullname = '%s [nMg=%d]' % (name, nMg_above)
                    #self.Add(groupvec, ddG0, fullname, cid, 'Mg')
                    # TODO: solve the problems with Mg species !!!
                    # there seems to be a problem in decomposition of ADP:Mg (and maybe other forms that have Mg ions)

    def AddFormationEnergies(self, obs_fname="../data/thermodynamics/formation_energies.csv"):
        dissociation = DissociationConstants.FromFile()
        train_species = PsuedoisomerTableThermodynamics.FromCsvFile(obs_fname, label='train')
        for cid in train_species.get_all_cids():
            pmap = train_species.cid2PseudoisomerMap(cid)
            pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
            if len(pmatrix) != 1:
                raise Exception("C%05d has more than one species in the training set" % cid)
            nH, charge, nMg, dG0 = pmatrix[0]
            name = "%s (%d)" % (self.kegg.cid2name(cid), nH)
            logging.info('Adding the formation energy of %s', name)
            self.html_writer.write("<h3>%s, %s</h3>\n" % (name, train_species.cid2SourceString(cid)))
            self.html_writer.write('&#x394;G<sub>f</sub> = %.2f, nH = %d, nMg = %d</br>\n' % (dG0, nH, nMg))
            smiles = dissociation.GetSmiles(cid, nH=nH, nMg=nMg)
            if not smiles:
                raise Exception("%s does not have a SMILES expression in the dissociation constant table" 
                                % name)
            self.html_writer.write("SMILES = %s</br>\n" % smiles)
            mol = Molecule.FromSmiles(smiles)
            mol.SetTitle(name)
            try:
                self.html_writer.write(mol.ToSVG())
            except (TypeError, IndexError, AssertionError):
                logging.warning('PyBel cannot draw the compound %s',  name)
                self.html_writer.write('WARNING: cannot draw this compound using PyBel\n')
    
            self.html_writer.write('</br>\n')
    
            try:
                decomposition = self.group_decomposer.Decompose(mol, strict=True)
            except GroupDecompositionError as e:
                logging.error('Cannot decompose one of the compounds in the training set: ' + str(mol))
                raise e
            
            groupvec = decomposition.AsVector()
            self.Add(groupvec, dG0, name, id="C%05d" % cid, obs_type='formation')
            self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
            
            gc_nH, gc_charge = decomposition.Hydrogens(), decomposition.NetCharge()
            if nH != gc_nH:
                s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, formula = %d' % (
                    nH, gc_nH)
                logging.error(s)
                self.html_writer.write(s + '</br>\n')
            if charge != gc_charge:
                s = 'ERROR: Charge doesn\'t match: explicit = %d, formula = %d' % (
                    charge, gc_charge)
                logging.error(s)
                self.html_writer.write(s + '</br>\n')
            
    def AddNistDatabase(self):
        nist_regression = NistRegression(self.db, html_writer=NullHtmlWriter())
        nist_regression.nist.override_pMg = 14.0
        cid2dG0 = {}
        cid2nH = {} # the nH that is to be used in the reverse transform
        for cid in nist_regression.dissociation.GetAllCids():
            diss = nist_regression.dissociation.GetDissociationTable(cid)
            nH, _nMg = diss.GetMostAbundantPseudoisomer(pH=default_pH, I=default_I, pMg=14, T=default_T)
            cid2nH[cid] = nH
        
        # override the nH in cid2nH with the pseudoisomer used in the CSV
        # so it will be consistent with the formation energy table
        train_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='train')
        for cid in train_species.get_all_cids():
            pmap = train_species.cid2PseudoisomerMap(cid)
            pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
            if len(pmatrix) != 1:
                raise Exception("C%05d has more than one species in the training set" % cid)
            cid2nH[cid] = pmatrix[0][0]

        # gather the dG0 for test compounds, so that their contribution
        # will be subtracted from the reaction dG0 and not used in the regression
        # of the groups
        test_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='test')
        for cid in test_species.get_all_cids():
            pmap = test_species.cid2PseudoisomerMap(cid)
            pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
            if len(pmatrix) != 1:
                raise Exception("C%05d has more than one species in the training set" % cid)
            nH, _z, _nMg, dG0 = pmatrix[0]
            cid2nH[cid] = nH
            cid2dG0[cid] = dG0

        S, dG0, cids = nist_regression.ReverseTransform(cid2nH=cid2nH)
        group_matrix = []
        good_cids = []
        good_indices = []
        anchored_cids = []
        dG0_correction = np.zeros((len(cids), 1))
        self.html_writer.write("<h3>The compounds used in NIST and their decompositions:</h3>\n")
        for i, cid in enumerate(cids):
            name = self.kegg.cid2name(cid)
            self.html_writer.write("<b>C%05d - %s</b></br>\n" % (cid, name))
            self.html_writer.write("nH = %d</br>\n" % cid2nH[cid])
            if cid in cid2dG0:
                # If the compound is marked for "test", normalize its contribution
                # and don't try to use it in the estimation
                self.html_writer.write("Anchored: dG0 = %.1f</br>\n" % 
                                       (cid2dG0[cid]))
                dG0_correction[i, 0] = cid2dG0[cid]
                anchored_cids.append(cid)
            else:
                try:
                    smiles = nist_regression.dissociation.GetSmiles(cid, nH=cid2nH[cid], nMg=0)
                    if not smiles:
                        raise GroupDecompositionError("no SMILES in the dissociation table")
        
                    self.html_writer.write("SMILES = %s</br>\n" % smiles)
                    mol = Molecule.FromSmiles(smiles)
                    self.html_writer.write(mol.ToSVG() + '</br>\n')
                
                    decomposition = self.group_decomposer.Decompose(mol, 
                        ignore_protonations=False, strict=True)
                    groupvec = decomposition.AsVector()
                    self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                    group_matrix.append(groupvec)
                    good_indices.append(i)
                    good_cids.append(cid)
                    gc_nH = decomposition.Hydrogens()
                    if cid2nH[cid] != gc_nH:
                        s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, formula = %d' % (
                            cid2nH[cid], gc_nH)
                        logging.error(s)
                        self.html_writer.write(s + '</br>\n')
                except GroupDecompositionError as e:
                    self.html_writer.write("ERROR:" + str(e) + "</br>\n")
        
        # The reverse transform converts the dG0' to dG0, where the chosen pseudoisomer
        # is the one with the minimal number of hydrogens. In order to correct this 
        # we need to subtract the difference between the minimal nH pseudoisomer and
        # the one used in the groupvector.
        dG0_noanchors = dG0 - np.dot(S, dG0_correction)
        
        # calculate the group vectors for each reaction, by multiplying the stoichiometric
        # matrix by the group matrix. Note that some rows (CIDs) are missing from the group_matrix
        # and the corresponding rows in S must also be skipped.
        observed_group_matrix = np.dot(S[:, good_indices], np.array(group_matrix))

        self.html_writer.write("<h3>The reaction used in NIST and their decompositions:</h3>\n")
        for r in xrange(S.shape[0]):
            name = "NIST%03d" % r
            sparse = dict([(cids[i], S[r, i]) for i in S[r, :].nonzero()[0]])

            self.html_writer.write("<p>\n<b>%s</b></br>\n" % name)
            self.html_writer.write("Original reaction: " + self.kegg.sparse_to_hypertext(sparse, show_cids=False) + "</br>\n")
            self.html_writer.write('&#x394;G<sub>r</sub> = %.2f</br>\n' % dG0[r, 0])

            for cid in set(sparse.keys()).intersection(anchored_cids):
                del sparse[cid]
            if sparse:
                min_cid = min(sparse.keys())
                if sparse[min_cid] > 0: # normalize the reaction direction such that the minimal CID will be on the left
                    sparse = dict([(cid, -coeff) for (cid, coeff) in sparse.iteritems()])
                    dG0_noanchors[r, 0] = -dG0_noanchors[r, 0]
                    observed_group_matrix[r, :] = -observed_group_matrix[r, :]
            
            self.html_writer.write("Truncated reaction: " + self.kegg.sparse_to_hypertext(sparse, show_cids=False) + "</br>\n")
            self.html_writer.write('&#x394;G<sub>r</sub> (truncated) = %.2f</br>\n' % dG0_noanchors[r, 0])

            unresolved_cids = set(sparse.keys()).difference(good_cids)
            if not unresolved_cids:
                groupvec = GroupVector(self.groups_data, observed_group_matrix[r, :])
                if groupvec: # if nonzero
                    self.Add(groupvec, dG0_noanchors[r, 0], name=name, id="", obs_type="reaction")
                    self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                else:
                    self.html_writer.write("Decompositions are equal on both sides - not using as an example</br>\n")
            else:
                self.html_writer.write("Some CIDs cannot be decomposed (%s) - not using as an example</br>\n" %
                                       ', '.join(['C%05d' % cid for cid in unresolved_cids]))
            self.html_writer.write("</p>\n")
        
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
                 ['dG0 REAL', 'name TEXT', 'id TEXT', 'obs_type TEXT']
        self.db.CreateTable('group_observations', ','.join(titles))

        self.db.CreateTable('pseudoisomer_observation', 
            'cid INT, name TEXT, protons INT, charge INT'
            ', nMg INT, dG0_f REAL, use_for TEXT')
        
        for obs in self.observations:
            obs.ToDatabase(self.db)
        self.db.Commit()
            
    def FromDatabase(self):
        self.observations = []
        for row in self.db.Execute("SELECT * FROM group_observations"):
            obs = GroupObservation()
            obs.groupvec = GroupVector(self.groups_data, row[0:self.n_groups])
            obs.dG0, obs.name, obs.cid, obs.obs_type = row[self.n_groups:]
            self.observations.append(obs)
        
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
            hash2obs_vec.setdefault(obs, []).append(obs)
        
        A = []
        b = []
        obs_types = []
        names = []
        for h, obs_vec in hash2obs_vec.iteritems():
            A.append(h.groupvec)
            b.append(np.mean([obs.dG0 for obs in obs_vec]))
            obs_types.append(h.obs_type)
            names.append(';'.join([obs.name for obs in obs_vec]))
            
        return np.matrix(A), np.array(b), obs_types, names
