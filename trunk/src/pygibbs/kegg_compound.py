#!/usr/bin/python

import logging
import pybel
import re

from pygibbs import elements
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
        self.num_electrons = None
        self.mol = None 
        self.from_kegg = True
        self.pubchem_id = None
        self.cas = ""

    def GetMolecule(self):
        if self.mol:
            return self.mol
        elif self.inchi:
            self.mol = Molecule.FromInChI(self.inchi)
            self.mol.SetTitle(self.name)
            return self.mol
        else:
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
        if row_dict['num_electrons']:
            comp.num_electrons = int(row_dict['num_electrons'])
        if row_dict['inchi']:
            comp.inchi = str(row_dict['inchi'])
            #comp.mol = Molecule.FromInChI(comp.inchi)
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
        return self.GetMolecule().ToInChI()
    
    def get_smiles(self):
        """Get a SMILES expression for this compound."""
        return self.GetMolecule().ToSmiles()

    def get_nH_and_charge(self):
        if not self.mol and self.inchi:
            self.mol = Molecule.FromInChI(self.inchi)

        if self.mol:
            return self.mol.GetHydrogensAndCharge()

        # if there is no InChI assume that the formula is correct and that
        # it represents the number of H for the neutral species
        atom_bag = self.get_atom_bag()
        if atom_bag:
            return atom_bag.get('H', 0), 0

        return None
        
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
        if self.mol:
            d['InChI'] = self.get_inchi()
            d['num_electrons'] = self.get_num_electrons()
        if self.mass is not None:
            d['mass'] = self.mass
        if self.formula:
            d['formula'] = self.formula
        if self.all_names:
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
