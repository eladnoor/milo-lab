#!/usr/bin/python

import logging
import re

from pygibbs import kegg_errors
from pygibbs import kegg_utils
from toolbox.molecule import Molecule

class Compound(object):
    """A class representing a compound in KEGG."""
    
    free_cid = -1 # class static variable
    def __init__(self, cid=None):
        """Initialize the Compound.
        
        Args:
            cid: the KEGG ID as an integer.
        
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
        self.name = "?"
        self.all_names = []
        self.mass = None
        self.formula = None
        self.inchi = None
        self.mol = None 
        self.from_kegg = True
        self.pubchem_id = None
        self.cas = ""

    def GetMolecule(self):
        """Gets a Molecule for this compound if possible.
        
        Returns None if no molecular data is available.
        """
        if self.mol:
            return self.mol
        
        if self.inchi:
            self.mol = Molecule.FromInChI(self.inchi)
            self.mol.SetTitle(self.name)
            return self.mol
        
        raise kegg_errors.KeggParseException(
             "C%05d doesn't have an explicit molecular structure" % self.cid)

    @staticmethod
    def FromDBRow(row_dict):
        """Build a Compound from a database row."""
        comp = Compound(cid=row_dict['cid'])
        comp.name = row_dict['name']
        comp.all_names = row_dict['all_names'].split(';')
        comp.mass = row_dict['mass']
        comp.formula = row_dict['formula']
        
        inchi = row_dict['inchi']
        if inchi is not None:
            comp.inchi = str(inchi)
        
        cas = row_dict['cas']
        if cas is not None:
            comp.cas = str(cas)
        
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
        if not atom_bag:
            return None
        
        atom_vector = [0] * Molecule.GetNumberOfElements()
        for elem, count in atom_bag.iteritems():
            if elem in ['R', 'X']:
                return None # wildcard compound!
            an = Molecule.GetAtomicNum(elem)
            if not an:
                logging.warning("Unsupported element in (C%05d): %s",
                                (self.cid, elem))
                return None
            atom_vector[an] = count
        return atom_vector
    
    def get_inchi(self):
        """Get the inchi code for this compound."""
        return self.GetMolecule().ToInChI()
    
    def get_smiles(self):
        """Get a SMILES expression for this compound."""
        return self.GetMolecule().ToSmiles()

    def get_nH_and_charge(self):
        if not self.mol and self.inchi:
            self.mol = Molecule.FromInChI(self.inchi)

        if self.mol:
            return self.mol.GetHydrogensAndCharge()

        # if there is no InChI assume that self.formula is correct and that
        # it represents the number of H for the neutral species
        atom_bag = self.get_atom_bag()
        if not atom_bag:
            return None
        return atom_bag.get('H', 0), 0
        
    def get_num_electrons(self):
        """Return the putative number of electrons in the molecule."""
        mol = self.GetMolecule()
        if mol:
            return mol.GetNumElectrons()
        
        # if there is no InChI assume that self.formula is correct and that
        # the charge is 0.
        atom_bag = self.get_atom_bag()
        if not atom_bag:
            return None
        n_protons = 0
        for elem, count in atom_bag.iteritems():
            n_protons += count * Molecule.GetAtomicNum(elem)
        return n_protons

    def get_link(self):
        """Returns a link to KEGG for this compound."""
        return kegg_utils.cid2link(self.cid)
    kegg_link = property(get_link)

    def ToJSONDict(self, verbose=False):
        """Converts to a JSON-formatted dictionary."""
        d = {}
        if self.cid:
            d['CID'] = "C%05d" % self.cid

        try:
            d['InChI'] = self.get_inchi()
        except Exception, e:
            d['InChI'] = None
        
        try:            
            d['num_electrons'] = self.get_num_electrons()
        except Exception, e:
            d['num_electrons'] = None
        
        if self.mass is not None:
            d['mass'] = self.mass
        if self.formula is not None:
            d['formula'] = self.formula
        if self.all_names is not None:
            d['names'] = self.all_names
        
        return d
    
    def __str__(self):
        return '%s: %s' % (self.cid, self.name)

def GetAllCompoundsFromDB(db):
    """Fetch all the compounds from the database."""
    compound_list = []
    for row_dict in db.DictReader('kegg_compound'):
        compound_list.append(Compound.FromDBRow(row_dict))
    return compound_list
