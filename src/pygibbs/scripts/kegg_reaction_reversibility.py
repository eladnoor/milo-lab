#!/usr/bin/python

import sys
from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics


def GetReactionIdInput():
    while True:
        try:
            print 'KEGG reaction ID:',
            return int(raw_input())
        except Exception:
            print 'KEGG reaction IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pH = options.ph
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg.getInstance()
    G = PsuedoisomerTableThermodynamics.FromDatabase(
        db, 'gc_pseudoisomers', name='milo_gc')
    
    print ('Parameters: T=%f K, pH=%.2g, pMg=%.2g, '
           'I=%.2gM, Median concentration=%.2gM' % (T, pH, pMg, I, c_mid))
    
    cmap = {}
    if not options.ignore_cofactors:
        print 'Fixing concentrations of co-factors'
        cmap = reversibility.GetConcentrationMap(kegg)
    else:
        print 'Not fixing concentrations of co-factors'

    while True:
        try:
            rid = GetReactionIdInput()        
            reaction = kegg.rid2reaction(rid)
            print 'Reaction Name: %s' % reaction.name
            print '\tKegg ID: R%05d' % rid
            print '\tEC: %s' % str(reaction.ec_list)
            rev = reversibility.CalculateReversability(reaction.sparse,
                                                       G, pH=pH, I=I, pMg=pMg,
                                                       T=T, concentration_map=cmap)
            print '\tIrreversibility: %.1f' % rev
        except Exception, e:
            print 'Error: ', e


if __name__ == '__main__':
    main()
    
