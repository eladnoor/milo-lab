import logging
import os
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.group_vector import GroupVector
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.nist_regression import NistRegression
from toolbox.molecule import Molecule
from toolbox.html_writer import HtmlWriter, NullHtmlWriter

class GroupObservation(object):
    
    def __init__(self):
        self.groupvec = None
        self.dG0 = None
        self.name = None
        self.id = None
        self.obs_type = None # can be 'formation', 'acid-base' or 'Mg'

    def NetCharge(self):
        return self.groupvec.NetCharge()

    def Hydrogens(self):
        return self.groupvec.Hydrogens()

    def Magnesiums(self):
        return self.groupvec.Magnesiums()

    def ToDatabase(self, db):
        db.Insert('group_observations', list(self.groupvec) +
                  [self.dG0, self.name, self.id, self.obs_type])
        
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
        
    def Add(self, groupvec, dG0, name, id, obs_type):
        obs = GroupObservation()
        obs.groupvec = groupvec
        obs.dG0 = dG0
        obs.name = name
        obs.id = id
        obs.obs_type = obs_type
        self.observations.append(obs)

    def SmilesTGroupvec(self, smiles, id):
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
            self.html_writer.write('%s Could not be decomposed</br>\n' % smiles)

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
            self.html_writer.write('%s: &#x394;G = %.2f</br>\n' 
                                   % (s, ddG0))
            self.html_writer.write('SMILES = %s >> %s</br>\n' % (smiles_below, smiles_above))
            decomposition_below = self.SmilesTGroupvec(smiles_below, 
                            "C%05d_b_H%d_Mg%d" % (cid, nH_below, nMg_below))
            decomposition_above = self.SmilesTGroupvec(smiles_above, 
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

    def AddPseudoisomersData(self, ps_isomer):
        #name = str(ps_isomer)
        name = "%s (%d)" % (ps_isomer.name, ps_isomer.net_charge)
        logging.info('Verifying data for %s', name)
        self.html_writer.write("<h3>%s, %s</h3>\n" % (ps_isomer.name, ps_isomer.ref))

        if ps_isomer.dG0 == None:
            self.html_writer.write('No data for &#x394;G<sub>f</sub></br>\n')
            return

        if ps_isomer.Skip():
            self.html_writer.write('Compound marked as not to be used</br>\n')
            return
            
        self.html_writer.write('&#x394;G<sub>f</sub> = %.2f</br>\n' % ps_isomer.dG0)

        if ps_isomer.cid:
            self.AddPseudoisomer(ps_isomer.cid, ps_isomer.hydrogens, ps_isomer.net_charge,
                                ps_isomer.magnesiums, ps_isomer.dG0)

        if ps_isomer.Test():
            self.cid_test_set.add(ps_isomer.cid)
            self.html_writer.write('Compound marked to be used only for testing (not training)</br>\n')
            return
        
        elif ps_isomer.Train():
            self.html_writer.write('Compound marked to be used for training</br>\n')
        else:
            raise Exception('Unknown usage flag: %' % ps_isomer.use_for)

        if not ps_isomer.smiles:
            raise Exception("Cannot use compound '%s' for training if it lacks a SMILES string" % ps_isomer.name)
        try:
            self.html_writer.write('SMILES = %s</br>\n' % ps_isomer.smiles)
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

        self.html_writer.write('</br>\n')

        try:
            decomposition = self.group_decomposer.Decompose(mol, strict=True)
        except GroupDecompositionError as e:
            logging.error('Cannot decompose one of the compounds in the training set: ' + str(mol))
            raise e
        
        groupvec = decomposition.AsVector()
        self.Add(groupvec, ps_isomer.dG0, name, id="C%05d" % ps_isomer.cid, obs_type='formation')
        self.html_writer.write("Decomposition = %s</br>\n" % decomposition)
        
        gc_hydrogens, gc_charge = decomposition.Hydrogens(), decomposition.NetCharge()
        if ps_isomer.hydrogens != gc_hydrogens:
            s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, formula = %d' % (
                ps_isomer.hydrogens, gc_hydrogens)
            logging.error(s)
            self.html_writer.write(s + '</br>\n')
        if ps_isomer.net_charge != gc_charge:
            s = 'ERROR: Charge doesn\'t match: explicit = %d, formula = %d' % (
                ps_isomer.net_charge, gc_charge)
            logging.error(s)
            self.html_writer.write(s + '</br>\n')
            
    def AddNistDatabase(self):
        nist_regression = NistRegression(self.db, html_writer=NullHtmlWriter())
        S, dG0, cids = nist_regression.ReverseTransform(use_anchors=False)
        
        group_matrix = []
        good_cids = []
        for cid in cids:
            name = self.kegg.cid2name(cid)
            mol = self.kegg.cid2mol(cid)
            min_nH, min_charge = self.kegg.cid2nH_and_charge(cid)
            self.html_writer.write("<h3>C%05d - %s</h3>\n" % (cid, name))
            self.html_writer.write("nH = %d, z = %d</br>\n" % (min_nH, min_charge))
            try:
                # TODO: the reverse transform converts the dG0' to dG0, where the chosen pseudoisomer
                # is the one with the minimal number of hydrogens. We must check that the nH
                # of the decomposition matches that of the entire compound.
                decomposition = self.group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
                groupvec = decomposition.AsVector()
                self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                if decomposition.Hydrogens() != min_nH:
                    self.html_writer.write("<b>ERROR</b>: decomposition nH (%d) is wrong</br>\n" % decomposition.Hydrogens())
                group_matrix.append(groupvec)
                good_cids.append(cid)
            except GroupDecompositionError as e:
                self.html_writer.write(str(e) + "</br>\n")
        
        # calculate the group vectors for each reaction, by multiplying the stoichiometric
        # matrix by the group matrix. Note that some rows (CIDs) are missing from the group_matrix
        # and the corresponding rows in S must also be skipped.
        good_indices = [cids.index(cid) for cid in good_cids]
        observed_group_matrix = np.dot(S[:, good_indices], np.array(group_matrix))

        for r in xrange(S.shape[0]):
            name = "NIST%03d" % r
            sparse = dict([(cids[i], S[r, i]) for i in S[r, :].nonzero()[0]])

            self.html_writer.write("<h3>%s</h3>\n" % name)
            self.html_writer.write('&#x394;G<sub>r</sub> = %.2f</br>\n' % dG0[r, 0])
            self.html_writer.write("Reaction: " + self.kegg.sparse_to_hypertext(sparse, show_cids=False) + "</br>\n")

            if set(sparse.keys()).issubset(good_cids):
                groupvec = GroupVector(self.groups_data, observed_group_matrix[r, :])
                if groupvec: # if nonzero
                    self.Add(groupvec, dG0[r, 0], name=name, id="", obs_type="reaction")
                    self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                else:
                    self.html_writer.write("Decompositions are equal on both sides - not using as an example</br>\n")
            else:
                self.html_writer.write("Some CIDs cannot be decomposed - not using as an example</br>\n")
        
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
            
            all_cids = '; '.join([obs.id or "" for obs in obs_vec])
            all_names = '; '.join([obs.name or "" for obs in obs_vec])
            if obs_vec[0].obs_type == 'formation':
                names.append(all_names)
            else:
                names.append('pK:(%s) - [%s]' % (h, all_cids))
            b.append(np.mean([obs.dG0 for obs in obs_vec]))
            
        return np.matrix(A), np.array(b), names
