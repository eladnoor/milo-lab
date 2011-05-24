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


def GetSparseReactionInput():
    while True:
        try:
            print 'KEGG reaction ID:',
            kegg_reaction = raw_input()        
            reaction = NistRowData.ParseReactionFormula(kegg_reaction, kegg_reaction)
            return reaction
        except KeggParseException:
            print 'KEGG reaction IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg.getInstance()
    thermo = PsuedoisomerTableThermodynamics.FromDatabase(
        db, 'gc_pseudoisomers', name='milo_gc')
    
    cmap = {}
    if not options.ignore_cofactors:
        print 'Fixing concentrations of co-factors'
        cmap = reversibility.GetConcentrationMap(kegg)
    else:
        print 'Not fixing concentrations of co-factors'

    print ('Parameters: T=%f K, pMg=%.2g, '
           'I=%.2gM, Median concentration=%.2gM' % (T, pMg, I, c_mid))
    
    while True:
        try:
            reaction = GetSparseReactionInput()
            print 'Reaction: %s' % reaction.HashableReactionString()
            pH_range = pylab.arange(5, 10.001, 0.1)
            dG0 = []
            for pH in pH_range:
                dG0.append(reaction.PredictReactionEnergy(thermo, pH=pH, pMg=pMg, I=I, T=T))
            pylab.plot(pH_range, dG0, '-')
            pylab.show()
                
        except Exception, e:
            print 'Error: ', e


if __name__ == '__main__':
    main()
    
