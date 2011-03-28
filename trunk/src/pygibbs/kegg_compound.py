#!/usr/bin/python

import logging
import openbabel
import pybel
import re

from pygibbs import elements
from pygibbs import kegg_errors
from pygibbs import kegg_utils


class Compound(object):
    """A class representing a compound in KEGG."""
    
    free_cid = -1 # class static variable
    
    def __init__(self, cid=None, name=None, all_names=None, mass=None,
                 formula=None, inchi=None):
        """Initialize the Compound.
        
        Args:
            cid: the KEGG ID as an integer.
            name: the primary name of the compound.
            all_names: the list of common names of the compound.
            mass: the molecular mass.
            formula: the chemical formula.
            inchi: the inchi code for the compound.
        
        Attributes:
            cid: the KEGG ID
            name: the primary name or "?"
            all_names: the list of common names all_names or []
            mass
            formula
            inchi
            from_kegg: whether the compound is from KEGG
            pubchem_id
            cas
        """
        if cid == None:
            self.cid = Compound.free_cid
            Compound.free_cid -= 1
        else:
            self.cid = cid
        self.name = name or "?"
        self.all_names = all_names or []
        self.mass = mass
        self.formula = formula
        self.num_electrons = None
        self.inchi = inchi
        self.from_kegg = True
        self.pubchem_id = None
        self.cas = ""
    
    @staticmethod
    def FromDBRow(row_dict):
        """Build a Compound from a database row."""
        comp = Compound(cid=row_dict['cid'])
        comp.name = row_dict['name']
        comp.all_names = row_dict['all_names'].split(';')
        comp.mass = row_dict['mass']
        comp.formula = row_dict['formula']
        comp.num_electrons = row_dict['num_electrons']
        if row_dict['inchi']:
            comp.inchi = str(row_dict['inchi'])
        if row_dict['cas']:
            comp.cas = str(row_dict['cas'])
        
        return comp
    
    def get_atom_bag(self):
        """Returns a dict containing the count for
           each atom in the compound.
        """
        if (self.formula == None or
            self.formula.find("(") != -1 or
            self.formula.find(")") != -1):
            return None
        
        atom_bag = {}
        for atom, count in re.findall("([A-Z][a-z]*)([0-9]*)", self.formula):
            if count == '':
                count = 1
            else:
                count = int(count)
            atom_bag[atom] = count

        if "R" in atom_bag:
            # This formula is not full ('R' is a wildcard not an atom)
            return None
        return atom_bag
    
    def get_atom_vector(self):
        atom_bag = self.get_atom_bag()
        atom_vector = [0] * len(elements.ELEMENTS.symbols)
        for elem, count in atom_bag.iteritems():
            if elem in ['R', 'X']:
                return None # wildcard compound!
            try:
                an = elements.ELEMENTS.symbol_to_an[elem]
                atom_vector[an-1] = count
            except KeyError:
                logging.warning(
                    "Unsupported element in (C%05d): %s", (self.cid, elem))
                return None
        return atom_vector
    
    def get_inchi(self):
        """Get the inchi code for this compound."""
        if not self.inchi:
            raise kegg_errors.KeggParseException(
                "C%05d doesn't have an 'inchi', so it is "
                "impossible get its molecular structure" % self.cid)
        return str(self.inchi)
    
    def get_smiles(self):
        """Get a SMILES expression for this compound."""
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("inchi", "smi")
        obmol = openbabel.OBMol()
        if not self.inchi:
            raise kegg_errors.KeggParseException(
                "C%05d doesn't have an 'inchi', so it "
                "cannot be converted it to SMILES" % self.cid)
        obConversion.ReadString(obmol, str(self.inchi))
        smiles_list = obConversion.WriteString(obmol).split()
        if smiles_list:
            return smiles_list[0]
        else:
            raise kegg_errors.KeggParseException(
                "C%05d cannot be convert to SMILES, because "
                "the InChI is %s" % str(self.inchi))

    def get_mol(self, remove_hydrogens=True):
        """
            I don't remember why it was necessary to first convert the string to SMILES and then
            create the Molecule object, but there must have been some kind of reason for it.
        """
        smiles = self.get_smiles()
        try:
            mol = pybel.readstring('smiles', smiles)
            if remove_hydrogens:
                mol.removeh()
        except IOError:
            raise kegg_errors.KeggParseException(
                "Cannot interpret the SMILES string for compound C%05d: %s"
                % (self.cid, smiles))
        mol.title = str(self.name)
        return mol
    
    def get_obmol(self, correctForPH=True, pH=7.4):
        if self.inchi == None:
            raise kegg_errors.KeggParseException(
                "C%05d doesn't have an 'inchi', "
                "so I can't get its molecular structure" % self.cid)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("inchi", "mol")
        obmol = openbabel.OBMol()
        obConversion.ReadString(obmol, self.get_inchi())
        if obmol.NumAtoms() == 0:
            raise kegg_errors.KeggParseException(
                "Cannot interpret the inchi string for compound C%05d: %s"
                % (self.cid, self.inchi))
        polaronly = False
        obmol.AddHydrogens(polaronly, correctForPH, pH)
        return obmol

    def get_nH(self, correctForPH=True, pH=7.4):
        """
            Returns the number of hydrogen atoms in a compound.
            It is calculated by subtracting the number of heavy atoms (anything bigger than H)
            from the total number of atoms.
        """
        try:
            obmol = self.get_obmol(correctForPH, pH)
            # HvyAtoms are all the non-hydrogen atoms
            return obmol.NumAtoms() - obmol.NumHvyAtoms()
        except kegg_errors.KeggParseException:
            atom_bag = self.get_atom_bag()
            if (atom_bag == None):
                return None
            else: 
                return atom_bag.get('H')
        
    def get_charge(self, correctForPH=True, pH=7.4):
        """Get the charge in the given pH."""
        try:
            return self.get_obmol(correctForPH, pH).GetTotalCharge()
        except kegg_errors.KeggParseException:
            return 0
    
    @staticmethod
    def CalculateNumElectrons(obmol):
        """Calculates the number of electrons in a given molecule."""
        atom_bag = {}
        for i in xrange(obmol.NumAtoms()):
            atom = obmol.GetAtom(i+1)
            atom_bag.setdefault(atom.GetAtomicNum(), 0)
            atom_bag[atom.GetAtomicNum()] += 1
        n_protons = sum([an*cnt for (an, cnt) in atom_bag.iteritems()])
        return n_protons - obmol.GetTotalCharge()
    
    def get_num_electrons(self):
        return self.num_electrons

    def get_link(self):
        """Returns a link to KEGG for this compound."""
        return kegg_utils.cid2link(self.cid)
    kegg_link = property(get_link)

    def ToJSONDict(self, verbose=False):
        d = {}
        if self.cid:
            d['CID'] = "C%05d" % self.cid
        if self.inchi:
            d['InChI'] = self.inchi
        if self.mass is not None:
            d['mass'] = self.mass
        if self.formula:
            d['formula'] = self.formula
        if self.all_names:
            d['names'] = self.all_names
        
        try:
            n_electrons = self.get_num_electrons()
            d['num_electrons'] = n_electrons
        except Exception, e:
            if verbose:
                logging.warning(e)
        
        return d
    
    def __str__(self):
        return '%s: %s' % (self.cid, self.name)

def GetAllCompoundsFromDB(db):
    """Fetch all the compounds from the database."""
    compound_list = []
    for row_dict in db.DictReader('kegg_compound'):
        compound_list.append(Compound.FromDBRow(row_dict))
    return compound_list
