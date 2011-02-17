#!/usr/bin/python

import logging
import pylab
import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import R, default_I, default_pH
from pygibbs.thermodynamic_constants import default_pMg, default_T
from pygibbs.groups import GroupContribution
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter


def GetReactionIdInput():
    while True:
        try:
            print 'KEGG reaction ID:',
            return int(raw_input())
        except Exception, e:
            print 'KEGG reaction IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pH = options.ph
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    G = GroupContribution(db, kegg=kegg)
    G.init()
    
    print ('Parameters: T=%f K, pH=%.2g, pMg=%.2g, '
           'I=%.2gmM, Median concentration=%.2gM' % (T, pH, pMg, I, c_mid))
    
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
            print 'Reaction Name', reaction.name
            print '\tKegg Id', reaction.rid
            print '\tEC', reaction.ec_list
            rev = reversibility.CalculateReversability(rid, G, pH=pH, I=I, pMg=pMg,
                                                       T=T, concentration_map=cmap)
            
            print '\tIrreversibility:', rev
        except Exception, e:
            print 'Error: ', e


if __name__ == '__main__':
    main()
    
