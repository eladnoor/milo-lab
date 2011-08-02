#!/usr/bin/python

import pylab
import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T, R
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from pygibbs.nist_verify import LoadAllEstimators

def GetReactionIdInput():
    while True:
        try:
            print 'KEGG compound ID:',
            return int(raw_input())
        except Exception:
            print 'KEGG compound IDs should be integers.'


def main():
    opt_parser = flags.MakeOpts()
    options, _ = opt_parser.parse_args(sys.argv)
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
    
    kegg = Kegg.getInstance()
    while True:
        cid = GetReactionIdInput()        
        compound = kegg.cid2compound(cid)
        print 'Compound Name: %s' % compound.name
        print '\tKegg ID: C%05d' % cid
        print '\tFormula: %s' % compound.formula
        print '\tInChI: %s' % compound.inchi
        for key, thermo in estimators.iteritems():
            print "\t<< %s >>" % key
            try:
                print thermo.cid2PseudoisomerMap(cid),
                print '--> dG0\'f = %.1f kJ/mol' % compound.PredictFormationEnergy(thermo)
            except Exception as e: 
                print '\t\tError: %s' % (str(e))


if __name__ == '__main__':
    main()
    
