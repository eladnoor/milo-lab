#!/usr/bin/python

import sys
from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.nist_verify import LoadAllEstimators


def GetReactionIdInput():
    while True:
        try:
            print 'KEGG reaction ID:',
            return int(raw_input())
        except Exception:
            print 'KEGG reaction IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    kegg = Kegg.getInstance()
    estimators = LoadAllEstimators()
    
    print ('Parameters: T=%f K, pH=%.2g, pMg=%.2g, '
           'I=%.2gmM, Median concentration=%.2gM' % 
           (default_T, options.ph, options.pmg, options.i_s, options.c_mid))

    for thermo in estimators.values():
        thermo.c_mid = options.c_mid
        thermo.pH = options.ph
        thermo.pMg = options.pmg
        thermo.I = options.i_s
        thermo.T = default_T
    
    cmap = {}
    if not options.ignore_cofactors:
        print 'Fixing concentrations of co-factors'
        cmap = reversibility.GetConcentrationMap()
    else:
        print 'Not fixing concentrations of co-factors'

    while True:
        rid = GetReactionIdInput()        
        reaction = kegg.rid2reaction(rid)
        print 'Reaction Name: %s' % reaction.name
        print '\tKegg ID: R%05d' % rid
        print '\tEC: %s' % str(reaction.ec_list)
        for key, thermo in estimators.iteritems():
            print "\t<< %s >>" % key
            try:
                print '\t\tdG0\'f = %.1f kJ/mol' % reaction.PredictReactionEnergy(thermo)
                rev = reversibility.CalculateReversability(reaction, thermo, concentration_map=cmap)
                print '\t\tgamma = %.3g' % rev
            except Exception as e: 
                print '\tError: %s' % (str(e))


if __name__ == '__main__':
    main()
    
