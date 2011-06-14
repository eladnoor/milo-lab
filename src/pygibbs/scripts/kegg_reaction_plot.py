#!/usr/bin/python

import sys
from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.nist import NistRowData
import pylab
from pygibbs.kegg_errors import KeggParseException
from pygibbs.hatzimanikatis import Hatzi


def GetSparseReactionInput():
    while True:
        try:
            print 'KEGG reaction ID:',
            kegg_reaction = raw_input()        
            reaction = NistRowData.ParseReactionFormula(kegg_reaction, kegg_reaction)
            return reaction
        except KeggParseException:
            print 'KEGG compound IDs should be integers (e.g. 2 + 1 = 8 + 9)'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase('../res/gibbs.sqlite')
    thermo_list = []
    thermo_list += [PsuedoisomerTableThermodynamics.FromDatabase(
        db, 'gc_pseudoisomers', name='milo_gc')]
    thermo_list += [Hatzi(use_pKa=False)]
    
    print ('Parameters: T=%f K, pMg=%.2g, '
           'I=%.2gM, Median concentration=%.2gM' % (T, pMg, I, c_mid))
    
    while True:
        try:
            reaction = GetSparseReactionInput()
            print 'Reaction: %s' % reaction.HashableReactionString()
            pH_range = pylab.arange(4, 10.001, 0.1)
            dG0_matrix = []
            for thermo in thermo_list:
                dG0 = []
                for pH in pH_range:
                    dG0.append(reaction.PredictReactionEnergy(thermo, pH=pH, pMg=pMg, I=I, T=T))
                dG0_matrix.append(dG0)
            dG0_matrix = pylab.array(dG0_matrix)
            pylab.plot(pH_range, dG0_matrix.T, '-')
            pylab.show()
                
        except Exception, e:
            print 'Error: ', e


if __name__ == '__main__':
    main()
    
