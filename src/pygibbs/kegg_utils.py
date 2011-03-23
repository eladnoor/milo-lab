#!/usr/bin/python

import openbabel
import re

from pygibbs import kegg_errors, kegg_reaction
from pylab import find
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T

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

def write_kegg_pathway(html_writer, reactions, fluxes):

    def write_reaction(prefix, reaction, flux=1):
        if (flux == 1):
            html_writer.write('%sR%05d&nbsp;&nbsp;%s<br>\n' % (prefix, reaction.rid, str(reaction)))
        else:
            html_writer.write('%sR%05d&nbsp;&nbsp;%s (x%g)<br>\n' % (prefix, reaction.rid, str(reaction), flux))
    
    html_writer.write('<p style="font-family: courier; font-size:10pt">')
    html_writer.write('ENTRY' + '&nbsp;'*7 + 'M-PATHOLOGIC<br>\n')
    html_writer.write('SKIP' + '&nbsp;'*8 + 'FALSE<br>\n')
    html_writer.write('NAME' + '&nbsp;'*8 + 'M-PATHOLOGIC<br>\n')
    html_writer.write('TYPE' + '&nbsp;'*8 + 'MARGIN<br>\n')
    html_writer.write('CONDITIONS' + '&nbsp;'*2 + 'pH=%g,I=%g,T=%g<br>\n' % 
                      (default_pH, default_I, default_T))
    html_writer.write('C_MID' + '&nbsp;'*7 + '0.0001<br>\n')
    for r in range(len(reactions)):
        if (r == 0):
            write_reaction('REACTION' + '&nbsp;'*4, reactions[r], fluxes[r])
        else:
            write_reaction('&nbsp;'*12, reactions[r], fluxes[r])
    html_writer.write('///<br></p>\n')
    
def write_module_to_html(html_writer, S, rids, fluxes, cids):
    reactions = []
    for r in xrange(S.shape[0]):
        sparse = dict([(cids[c], S[r,c]) for c in find(S[r,:])])
        reaction = kegg_reaction.Reaction('R%05d' % rids[r], sparse, rid=rids[r])
        reactions.append(reaction)
    write_kegg_pathway(html_writer, reactions, fluxes)