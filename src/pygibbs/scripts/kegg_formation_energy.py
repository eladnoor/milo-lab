#!/usr/bin/python

import pylab
import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T, R
from pygibbs import reversibility
from pygibbs.kegg import Kegg
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase
from pygibbs.nist_verify import LoadAllEstimators

def GetReactionIdInput():
    while True:
        try:
            print 'KEGG compound ID:',
            return int(raw_input())
        except Exception:
            print 'KEGG compound IDs should be integers.'


def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    c_mid = options.c_mid
    pH = options.ph
    pMg = options.pmg
    I = options.i_s
    T = default_T

    db = SqliteDatabase("../res/gibbs.sqlite")
    observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
        '../data/thermodynamics/dG0.csv')
    gc_table_name = "gc_pseudoisomers"
    if not db.DoesTableExist(gc_table_name):
        raise ValueError('The table %s does not exist in the database. '
                         'Please run the groups.py script and try again.'
                         % gc_table_name)
    thermo = PsuedoisomerTableThermodynamics.FromDatabase(
        db, gc_table_name)
    thermo.override_data(observed_thermo)
    kegg = Kegg.getInstance()
    
    estimators = LoadAllEstimators()
    estimators['old estimate'] = thermo
    
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
            cid = GetReactionIdInput()        
            compound = kegg.cid2compound(cid)
            print 'Compound Name: %s' % compound.name
            print '\tKegg ID: C%05d' % cid
            print '\tFormula: %s' % compound.formula
            print '\tInChI: %s' % compound.inchi
            conc = cmap.get(cid, c_mid)
            print '\tConcentration: %.2e' % conc
            for key, thermo in estimators.iteritems():
                dG0_tag = thermo.cid2PseudoisomerMap(cid).Transform(pH, pMg, I, T)
                dG_tag = dG0_tag + R*T*pylab.log(conc)
                print '\t%s estimated dG0\'f at 1M (and at %gM): %.1f (%.1f)' % (key, conc, dG0_tag, dG_tag)
        except Exception, e:
            print 'Error: ', e


if __name__ == '__main__':
    main()
    
