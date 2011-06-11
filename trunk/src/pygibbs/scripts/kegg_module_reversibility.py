#!/usr/bin/python

import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs.groups import GroupContribution
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase


def GetModuleIdInput():
    while True:
        try:
            print 'KEGG module ID:',
            return int(raw_input())
        except Exception:
            print 'KEGG module IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pH = options.ph
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg.getInstance()
    G = GroupContribution(db)
    G.init()

    print ('Parameters: T=%f K, pH=%.2g, pMg=%.2g, '
           'I=%.2gM, Median concentration=%.2gM' % (T, pH, pMg, I, c_mid))
    
    cmap = {}
    if not options.ignore_cofactors:
        if options.full_metabolites:
            print 'Fixing concentrations of all known metabolites'
            cmap = reversibility.GetFullConcentrationMap(G)
        else:
            print 'Fixing concentrations of co-factors'
            cmap = reversibility.GetConcentrationMap(kegg)
    else:
        print 'Not fixing concentrations of co-factors'

    if options.report_mode:
        print 'Output used metabolites concentrations'

    while True:
        mid = GetModuleIdInput()
            
        rid_flux_list = kegg.mid2rid_map[mid]

        for rid, flux in rid_flux_list:
            try:
                reaction = kegg.rid2reaction(rid)
                print 'Reaction Name', reaction.name
                print '\tKegg Id', reaction.rid
                print '\tEC', reaction.ec_list
                rev = reversibility.CalculateReversability(reaction.sparse,
                                                           G, pH=pH, I=I, pMg=pMg,
                                                           T=T, concentration_map=cmap)
                if rev == None:
                    dG = G.estimate_dG_reaction(reaction.sparse, pH=pH, pMg=pMg, I=I, T=T, c0=c_mid, media='glucose')
                    print '\tReversibility: No free compounds, dG = %.2g' % dG
                else:
                    corrected_reversibility = flux * rev
                    print '\tReversibility %.2g' % corrected_reversibility
                    
                if options.report_mode:
                    for cid,s in reaction.sparse.iteritems():
                        if cid in cmap:
                            print '(%d C%05d) %s\t: %.2g' % (s, cid, kegg.cid2name(cid), cmap[cid])
                        else: 
                            print '(%d C%05d) %s\t: Free concentration' % (s, cid, kegg.cid2name(cid))
            except Exception:
                print '\tCouldn\'t calculate irreversibility'


if __name__ == '__main__':
    main()