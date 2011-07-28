#!/usr/bin/python

import pylab, logging
from pygibbs import kegg_errors
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T

##
# TODO: (flamholz): Not all these utilities are specific to KEGG.
##

def write_kegg_pathway(html_writer, reactions, fluxes):

    def write_reaction(prefix, reaction, flux=1):
        if (flux == 1):
            html_writer.write('%sR%05d&nbsp;&nbsp;%s<br>\n' % (prefix, reaction.rid,
                                                               reaction.FullReactionString()))
        else:
            html_writer.write('%sR%05d&nbsp;&nbsp;%s (x%g)<br>\n' % (prefix, reaction.rid,
                                                                     reaction.FullReactionString(), flux))
    
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
    from pygibbs.kegg_reaction import Reaction
    
    reactions = []
    for r in xrange(S.shape[0]):
        sparse = dict([(cids[c], S[r,c]) for c in pylab.find(S[r,:])])
        reaction = Reaction('R%05d' % rids[r], sparse, rid=rids[r])
        reactions.append(reaction)
    write_kegg_pathway(html_writer, reactions, fluxes)
    
def balance_reaction(kegg, sparse, balance_water=False, balance_hydrogens=False,
                     exception_if_unknown=False):
    """
        Balances a reaction
        
        Arguments:
            If balance_water=True and there is an imbalance of oxygen atoms, Balance
                changes the reaction by adding H2O until it is balanced.
            
            If balance_hydrogens=True then H+ are used to balance the amount of hydrogen atoms.
            
            If exception_if_unknown=True then an exception will be raised also if there
                is not enough information to know if the reaction is balanced or not.
        
        If the reaction cannot be balanced, raises KeggReactionNotBalancedException
    """
    atom_bag = {}
    try:
        for cid, coeff in sparse.iteritems():
            comp = kegg.cid2compound(cid)
            cid_atom_bag = comp.get_atom_bag()
            if cid_atom_bag == None:
                if exception_if_unknown:
                    raise kegg_errors.KeggReactionNotBalancedException(
                        "C%05d has no explicit formula, "
                        "cannot check if this reaction is balanced" % cid)
                else:
                    logging.warning("C%05d has no explicit formula, "
                                    "cannot check if this reaction is balanced" % cid)
                    return
            try:
                cid_atom_bag['e-'] = comp.get_num_electrons()
            except kegg_errors.KeggParseException:
                return
            
            for atomicnum, count in cid_atom_bag.iteritems():
                atom_bag[atomicnum] = atom_bag.get(atomicnum, 0) + count*coeff
                
    except KeyError as e:
        if exception_if_unknown:
            raise kegg_errors.KeggReactionNotBalancedException(
                "cannot check if this reaction is balanced")
        else:
            logging.warning(str(e) + ", cannot check if this reaction is balanced")
            return

    if balance_water and atom_bag.get('O', 0) != 0:
        sparse[1] = sparse.get(1, 0) - atom_bag['O'] # balance the number of oxygens by adding C00001 (water)
        atom_bag['H'] = atom_bag.get('H', 0) - 2 * atom_bag['O'] # account for the 2 hydrogens in each added water molecule
        atom_bag['e-'] = atom_bag.get('e-', 0) - 10 * atom_bag['O'] # account for the 10 electrons in each added water molecule
        atom_bag['O'] = 0
    
    if balance_hydrogens:
        if atom_bag.get('H', 0) != 0:
            sparse[80] = sparse.get(80, 0) - atom_bag['H'] # balance the number of hydrogens by adding C00080 (H+)
            atom_bag['H'] = 0
    else:
        if 80 in sparse:
            del sparse[80]
        atom_bag['H'] = 0
    
    for atomtype in atom_bag.keys():
        if atom_bag[atomtype] == 0:
            del atom_bag[atomtype]

    if atom_bag:
        raise kegg_errors.KeggReactionNotBalancedException("Reaction cannot be balanced: " 
            + str(atom_bag))
