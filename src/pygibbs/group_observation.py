import logging
import os
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.group_vector import GroupVector
from pygibbs.nist_regression import NistRegression
from toolbox.html_writer import NullHtmlWriter
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T,\
    default_pMg
import csv
from pygibbs.nist import Nist

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
        db.Insert('group_observations', self.ToCsvRow())
    
    @staticmethod 
    def FromDatabaseRow(row, groups_data):
        n_groups = len(groups_data.all_groups)
        obs = GroupObservation()
        obs.groupvec = GroupVector(groups_data, row[0:n_groups])
        obs.dG0, obs.name, obs.cid, obs.obs_type = row[n_groups:]
        return obs
        
    def ToCsvRow(self):
        return list(self.groupvec) + [self.dG0, self.name, self.id, self.obs_type]
    
    def __hash__(self):
        return hash(str(self.groupvec))
    
    def __eq__(self, other):
        return self.groupvec == other.groupvec
    
    def Normalize(self):
        if not self.groupvec.__nonzero__():
            return
        for v in self.groupvec:
            if v != 0:
                self.groupvec = self.groupvec * (1.0 / v)
                self.dG0 *= (1.0 / v)
                return
        
class GroupObervationCollection(object):
    
    def __init__(self, db, html_writer, group_decomposer, dissociation,
                 transformed=False):
        self.kegg = Kegg.getInstance()
        self.db = db
        self.groups_data = group_decomposer.groups_data
        self.group_decomposer = group_decomposer
        self.observations = []
        self.n_groups = len(self.groups_data.all_groups)
        self.groupvec_titles = ["g%d" % i for i in xrange(self.n_groups)]
        self.dissociation = dissociation
        self.transformed = transformed
        self.use_pkas = False
        self.html_writer = html_writer
        self.TRAIN_GROUPS_TABLE_NAME = 'pgc_train_groups'
        self.OBSERVATION_TABLE_NAME = 'pgc_observations'
        
    def FromFiles(self):
        self.html_writer.write('</br><b>List of NIST reactions for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        if not self.transformed:
            self.AddReverseTransformedNistDatabase()
        else:
            self.AddNistDatabase()
        self.html_writer.div_end()
        
        self.html_writer.write('<br><b>List of compounds for training</b>')
        self.html_writer.insert_toggle(start_here=True)
        if not self.transformed:
            self.AddChemicalFormationEnergies()
        else:
            self.AddBiochemicalFormationEnergies()
        self.html_writer.div_end()
        
        if self.use_pkas:
            self.html_writer.write('<b>List of pKa for training</b>')
            self.html_writer.insert_toggle(start_here=True)
            self.obs_collection.AddDissociationTable()
            self.html_writer.div_end()
        
    def Add(self, groupvec, dG0, name, id, obs_type):
        obs = GroupObservation()
        obs.groupvec = groupvec
        obs.dG0 = dG0
        obs.name = name
        obs.id = id
        obs.obs_type = obs_type
        obs.Normalize()
        self.observations.append(obs)

    def MolToGroupvec(self, mol, id):
        mol.RemoveHydrogens()
        mol.SetTitle(id)
        
        #try:
        self.html_writer.write(mol.ToSVG())
        #except (TypeError, AssertionError): # The Mg ions cannot be drawn by OASA 
        #    pass
        
        try:
            return self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
        except GroupDecompositionError as _e:
            logging.warning('Cannot decompose one of the compounds in the '
                            'training set: %s, %s' % (id, str(mol)))
            self.html_writer.write('%s Could not be decomposed</br>\n' % str(mol))

    def AddDissociationTable(self):
        for cid in self.dissociation.GetAllCids():
            diss_table = self.dissociation.GetDissociationTable(cid)
            if diss_table is None:
                continue
            name = self.kegg.cid2name(cid)
            logging.debug("Reading pKa data for C%05d, %s:" % (cid, name))
            self.html_writer.write('<h3>C%05d - %s</h3>\n' % (cid, name))
            for key in diss_table:
                nH_above, nH_below, nMg_above, nMg_below = key
                ddG0, _ref = diss_table.ddGs[key]
                mol_above = diss_table.GetMol(nH_above, nMg_above)
                mol_below = diss_table.GetMol(nH_below, nMg_below)
                
                if mol_above is None or mol_below is None:
                    continue
        
                if nH_above != nH_below:
                    s = "nH = %d -> %d" % (nH_below, nH_above)
                elif nMg_above != nMg_below:
                    s = "nMg = %d -> %d" % (nMg_below, nMg_above)
                logging.debug("\t" + s)
                self.html_writer.write('%s: &#x394;G = %.2f</br>\n' 
                                       % (s, ddG0))
                self.html_writer.write('SMILES = %s >> %s</br>\n' %
                                (mol_above.ToSmiles(), mol_below.ToSmiles()))
                decomposition_below = self.MolToGroupvec(mol_below, 
                                "C%05d_b_H%d_Mg%d" % (cid, nH_below, nMg_below))
                decomposition_above = self.MolToGroupvec(mol_above, 
                                "C%05d_a_H%d_Mg%d" % (cid, nH_above, nMg_above))
                if decomposition_below is None or decomposition_above is None:
                    continue
                groupvec = decomposition_below.AsVector() - decomposition_above.AsVector()
                self.html_writer.write('</br>\nDecomposition = %s</br>\n' % str(groupvec))
                
                if nH_above != nH_below:
                    if nH_above != nH_below-1:
                        raise Exception("a pKa must represent a difference of exactly 1 hydrogen")
                    fullname = '%s [nH=%d]' % (name, nH_above)
                    self.Add(groupvec, ddG0, fullname, "C%05d" % cid, obs_type='acid-base')
                elif nMg_above != nMg_below:
                    if nMg_above != nMg_below-1:
                        raise Exception("a pK_Mg must represent a difference of exactly 1 magnesium")
                    fullname = '%s [nMg=%d]' % (name, nMg_above)
                    #self.Add(groupvec, ddG0, fullname, cid, obs_type='Mg')
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
            
            name = "%s, %s" % (self.kegg.cid2name(cid), ref)
            logging.debug('Adding the formation energy of %s', name)
            self.html_writer.write("<b>%s (C%05d), %s</b></br>\n" % 
                                   (name, cid, ref))
            self.html_writer.write('&#x394;<sub>f</sub> G\'<super>o</super> = %.2f, '
                                   'pH = %g, I = %g, pMg = %g, T = %g</br>\n' % 
                                   (dG0_prime, pH, I, pMg, T))
            mol = self.dissociation.GetMostAbundantMol(cid, pH, I, pMg, T)
            if mol is None:
                raise Exception("%s [C%05d] does not have a SMILES "
                                "expression in the dissociation constant table" 
                                % (name, cid))
            self.html_writer.write('<font size="1">SMILES = %s</font></br>\n' % mol.ToSmiles())
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
            self.Add(groupvec, dG0_prime, name, id="C%05d" % cid, obs_type='formation')
            self.html_writer.write('<font size="1">Decomposition = %s</font></br>\n' % groupvec)
            
    def AddChemicalFormationEnergies(self,
            obs_fname="../data/thermodynamics/formation_energies.csv"):
        train_species = PsuedoisomerTableThermodynamics.FromCsvFile(obs_fname, label='training')
        for cid in train_species.get_all_cids():
            pmap = train_species.cid2PseudoisomerMap(cid)
            pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
            if len(pmatrix) != 1:
                raise Exception("multiple training species for C%05d" % cid)
            nH, charge, nMg, dG0 = pmatrix[0]
            name = "%s (%d)" % (self.kegg.cid2name(cid), nH)
            logging.debug('Adding the formation energy of %s', name)
            self.html_writer.write("<b>%s (C%05d), %s</b></br>\n" % 
                                   (name, cid, train_species.cid2SourceString(cid)))
            self.html_writer.write('&#x394;G<sub>f</sub> = %.2f, '
                                   'nH = %d, nMg = %d</br>\n' % (dG0, nH, nMg))
            diss_table = self.dissociation.GetDissociationTable(cid, 
                                                        create_if_missing=True)
            if diss_table is None:
                raise Exception("%s [C%05d, nH=%d, nMg=%d] does not have a " 
                                "dissociation table"
                                % (name, cid, nH, nMg))
            mol = diss_table.GetMol(nH=nH, nMg=nMg)
            if mol is None:
                raise Exception("%s [C%05d, nH=%d, nMg=%d] does not have a SMILES "
                                "expression in the dissociation constant table" 
                                % (name, cid, nH, nMg))
            self.html_writer.write('<font size="1">SMILES = %s</font></br>\n' % mol.ToSmiles())
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
            self.html_writer.write('<font size="1">Decomposition = %s</font></br>\n' % groupvec)
            
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
        nist = Nist()
        all_cids = sorted(nist.GetAllCids())

        S = np.zeros((0, len(all_cids))) # stoichiometric matrix
        dG0_prime = np.zeros((0, 1))
        for nist_row_data in nist.SelectRowsFromNist():
            dG0_prime = np.vstack([dG0_prime, nist_row_data.dG0_r])
            #pH = nist_row_data.pH
            #I = nist_row_data.I
            #pMg = nist_row_data.pMg
            stoichiometric_row = np.zeros((1, len(all_cids)))
            for cid, coeff in nist_row_data.reaction.iteritems():
                stoichiometric_row[0, all_cids.index(cid)] = coeff
            S = np.vstack([S, stoichiometric_row])

        # remove the columns that are all-zeros in S
        nonzero_columns = np.sum(abs(S), 0).nonzero()[0]
        cids = [all_cids[i] for i in nonzero_columns]
        S = S[:, nonzero_columns]

        cid2dG0_prime = {}
        obs_fname = '../data/thermodynamics/formation_energies_transformed.csv'
        for row in csv.DictReader(open(obs_fname, 'r')):
            cid = int(row['cid'])
            #pH = float(row['pH'])
            #I = float(row['I'])
            #pMg = float(row['pMg'])
            #T = float(row['T'])
            #ref = row['compound_ref']
            if row['use for'] == 'testing':
                cid2dG0_prime[cid] = float(row["dG'0"])

        group_matrix = []
        good_cids = []
        good_indices = []
        anchored_cids = []
        dG0_correction = np.zeros((len(cids), 1))
        self.html_writer.write("<h3>The compounds used in NIST and their decompositions:</h3>\n")
        for i, cid in enumerate(cids):
            name = self.kegg.cid2name(cid)
            self.html_writer.write("<b>C%05d - %s</b></br>\n" % (cid, name))
            if cid in cid2dG0_prime:
                # If the compound is marked for "test", normalize its contribution
                # and don't try to use it in the estimation
                self.html_writer.write("Anchored: dG0 = %.1f</br>\n" % 
                                       (cid2dG0_prime[cid]))
                dG0_correction[i, 0] = cid2dG0_prime[cid]
                anchored_cids.append(cid)
            else:
                try:
                    mol = self.dissociation.GetMostAbundantMol(cid, pH=default_pH, 
                        I=0.0, pMg=14.0, T=default_T)
                    if mol is None:
                        raise GroupDecompositionError("This pseudoisomer does not "
                        "have an explicit formula in the dissociation table: "
                        "C%05d" % cid)                

                    self.html_writer.write("SMILES = %s</br>\n" % mol.ToSmiles())
                    self.html_writer.write(mol.ToSVG() + '</br>\n')
                    decomposition = self.group_decomposer.Decompose(mol, 
                        ignore_protonations=False, strict=True)
                    groupvec = decomposition.AsVector()
                    self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                    group_matrix.append(groupvec)
                    good_indices.append(i)
                    good_cids.append(cid)
                except GroupDecompositionError as e:
                    self.html_writer.write("ERROR:" + str(e) + "</br>\n")
        
        dG0_prime_noanchors = dG0_prime - np.dot(S, dG0_correction)
        observed_group_matrix = np.dot(S[:, good_indices], np.array(group_matrix))
        self.html_writer.write("<h3>The reaction used in NIST and their decompositions:</h3>\n")
        for r in xrange(S.shape[0]):
            name = "NIST%03d" % r
            sparse = dict([(cids[i], S[r, i]) for i in S[r, :].nonzero()[0]])

            self.html_writer.write("<p>\n<b>%s</b></br>\n" % name)
            biochemical_reaction = " + ".join(['%g <a href="%s" title="C%05d">%s</a>' % 
                (coeff, self.kegg.cid2link(cid), cid, self.kegg.cid2name(cid)) 
                for (cid, coeff) in sparse.iteritems()])
            self.html_writer.write("Biochemical reaction: " + biochemical_reaction + "</br>\n")
            self.html_writer.write('&#x394;<sub>r</sub>G\'<super>o</super> = %.2f</br>\n' % dG0_prime[r, 0])

            for cid in set(sparse.keys()).intersection(anchored_cids):
                del sparse[cid]
            
            self.html_writer.write("Truncated reaction: " + self.kegg.sparse_to_hypertext(sparse, show_cids=False) + "</br>\n")
            self.html_writer.write('&#x394;<sub>r</sub>G\'<super>o</super> (truncated) = %.2f</br>\n' % dG0_prime_noanchors[r, 0])

            unresolved_cids = set(sparse.keys()).difference(good_cids)
            if not unresolved_cids:
                groupvec = GroupVector(self.groups_data, observed_group_matrix[r, :])
                if groupvec: # if nonzero
                    self.Add(groupvec, dG0_prime_noanchors[r, 0], name=name, id="", obs_type="reaction")
                    self.html_writer.write("Decomposition = %s</br>\n" % groupvec)
                else:
                    self.html_writer.write("Decompositions are equal on both sides - not using as an example</br>\n")
            else:
                self.html_writer.write("Some CIDs cannot be decomposed (%s) - not using as an example</br>\n" %
                                       ', '.join(['C%05d' % cid for cid in unresolved_cids]))
            self.html_writer.write("</p>\n")

    def AddReverseTransformedNistDatabase(self):
        nist_regression = NistRegression(self.db, dissociation=self.dissociation,
                                         html_writer=NullHtmlWriter())
        #nist_regression.nist.override_pMg = 14.0
        cid2dG0 = {}
        cid2nH = {} # the nH that is to be used in the reverse transform
        for cid in self.dissociation.GetAllCids():
            diss_table = self.dissociation.GetDissociationTable(cid)
            if diss_table is not None:
                nH, _nMg = diss_table.GetMostAbundantPseudoisomer(
                            pH=default_pH, I=0, pMg=14, T=default_T)
                cid2nH[cid] = nH
            else:
                # assume nH=0 by default for compounds without explicit formulas
                cid2nH[cid] = 0
        
        # override the nH in cid2nH with the pseudoisomer used in the CSV
        # so it will be consistent with the formation energy table
        train_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/formation_energies.csv', label='training')
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
            '../data/thermodynamics/formation_energies.csv', label='testing')
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
                    mol = self.dissociation.GetMol(cid, nH=cid2nH[cid], nMg=0)
                    if mol is None:
                        raise GroupDecompositionError("This pseudoisomer does not "
                        "have an explicit formula in the dissociation table: "
                        "C%05d [nH=%d, nMg=%d]" % (cid, cid2nH[cid], 0))
        
                    self.html_writer.write("SMILES = %s</br>\n" % mol.ToSmiles())
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
            chemical_reaction = " + ".join(['%g <a href="%s" title="C%05d">%s</a> [nH=%d]' % 
                (coeff, self.kegg.cid2link(cid), cid, self.kegg.cid2name(cid), cid2nH[cid]) 
                for (cid, coeff) in sparse.iteritems()])
            self.html_writer.write("Chemical reaction: " + chemical_reaction + "</br>\n")
            self.html_writer.write('&#x394;G<sub>r</sub> = %.2f</br>\n' % dG0[r, 0])

            for cid in set(sparse.keys()).intersection(anchored_cids):
                del sparse[cid]
            
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
        self.db.CreateTable(self.TRAIN_GROUPS_TABLE_NAME, 'name TEXT, nH INT, z INT, nMg INT')
        for group in self.groups_data.all_groups:
            self.db.Insert(self.TRAIN_GROUPS_TABLE_NAME, [group.name, group.hydrogens,
                                            group.charge, group.nMg])

        # the table 'group_observation' will contain all the observed data
        # that is used for training later
        titles = ['%s REAL' % t for t in self.groupvec_titles] + \
                 ['dG0 REAL', 'name TEXT', 'id TEXT', 'obs_type TEXT']
        self.db.CreateTable(self.OBSERVATION_TABLE_NAME, ','.join(titles))

        for obs in self.observations:
            obs.ToDatabase(self.db)
        self.db.Commit()
        
    def ToCSV(self, fname):
        csv_writer = csv.writer(open(fname, 'w'))
        titles = ['%s REAL' % t for t in self.groupvec_titles] + \
                 ['dG0 REAL', 'name TEXT', 'id TEXT', 'obs_type TEXT']
        csv_writer.writerow(titles)
        for obs in self.observations:
            csv_writer.writerow(obs.ToCsvRow())
            
    def FromDatabase(self):
        self.observations = []
        for row in self.db.Execute("SELECT * FROM %s" % self.OBSERVATION_TABLE_NAME):
            obs = GroupObservation.FromDatabaseRow(row, self.groups_data)
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
            names.append(' ; '.join([obs.name for obs in obs_vec]))
            
        return np.matrix(A), np.array(b), obs_types, names
