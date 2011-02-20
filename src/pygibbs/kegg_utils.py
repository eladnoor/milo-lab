#!/usr/bin/python

import openbabel
import re

from pygibbs import kegg_errors

##
## TODO(flamholz): Not all these utilities are specific to KEGG.
## Maybe move some of them elsewhere.
##


def mol2inchi(mol):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "inchi")
    return obConversion.WriteString(mol.OBMol).strip()


def mol2smiles(mol):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("inchi", "smi")
    return obConversion.WriteString(mol.OBMol).split()[0]


def smiles2inchi(smiles):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "inchi")
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, str(smiles))
    return obConversion.WriteString(obmol).strip()


def inchi2smiles(inchi):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("inchi", "smi")
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, inchi)
    return obConversion.WriteString(obmol).split()[0]


def remove_atoms_from_mol(mol, atoms):
    obmol = mol.OBMol
    obmol.BeginModify()
    for i in sorted(atoms, reverse=True):
        obmol.DeleteAtom(obmol.GetAtom(i))
    obmol.EndModify()


def cid2link(cid):
    """Returns the KEGG link for this compound."""
    return "http://www.genome.jp/dbget-bin/www_bget?cpd:C%05d" % cid


def parse_reaction_formula_side(s):
    """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        return the set of CIDs, ignore stoichiometry
    """
    if s.strip() == "null":
        return {}
    
    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            try:
                amount = float(tokens[0])
            except ValueError:
                raise kegg_errors.KeggParseException(
                    "Non-specific reaction: %s" % s)
            key = tokens[1]
            
        if key[0] != 'C':
            raise kegg_errors.KeggNonCompoundException(
                "Compound ID doesn't start with C: %s" % key)
        try:
            cid = int(key[1:])
            compound_bag[cid] = compound_bag.get(cid, 0) + amount
        except ValueError:
            raise kegg_errors.KeggParseException(
                "Non-specific reaction: %s" % s)
    
    return compound_bag


def parse_reaction_formula(formula):
    """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
        return the set of substrates, products and the direction of the reaction
    """
    tokens = re.findall("([^=^<]+) (<*=>*) ([^=^>]+)", formula)
    if len(tokens) != 1:
        raise kegg_errors.KeggParseException(
            "Cannot parse this formula: %s" % formula)
    
    left, direction, right = tokens[0] # the direction: <=, => or <=>
    
    sparse_reaction = {}
    for cid, count in parse_reaction_formula_side(left).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

    for cid, count in parse_reaction_formula_side(right).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

    return (sparse_reaction, direction)


def unparse_reaction_formula(sparse, direction='=>'):
    s_left = []
    s_right = []
    for cid, count in sparse.iteritems():
        show_string = "C%05d" % cid
        
        if count > 0:
            if count == 1:
                s_right.append(show_string)
            else:
                s_right.append('%d %s' % (count, show_string))
        elif count < 0:
            if count == -1:
                s_left.append(show_string)
            else:
                s_left.append('%d %s' % (-count, show_string))
    return ' + '.join(s_left) + ' ' + direction + ' ' + ' + '.join(s_right)