"""
    This script analyses the predictions of all the different estimation methods:
    Alberty, Hatzimanikatis, and the Milo lab Group Contribution method
"""

from pylab import * #@UnusedWildImport
from toolbox.html_writer import HtmlWriter
from pygibbs.groups import GroupContribution
from pygibbs.kegg import Kegg
from pygibbs.alberty import Alberty
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from toolbox import database
import logging
from pygibbs.nist_regression import NistRegression
from pygibbs.thermodynamics import CsvFileThermodynamics

################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter("../res/nist/report.html")
    kegg = Kegg()

    nist = Nist(db, html_writer, kegg)
    if True:
        nist.FromCsvFile('../data/thermodynamics/nist.csv')
    else:
        nist.FromDatabase()

    estimators = {}
    estimators['Alberty'] = CsvFileThermodynamics('../data/thermodynamics/alberty_pseudoisomers.csv')
    estimators['Hatzimanikatis Group Contribution'] = Hatzi()

    regress = NistRegression(db, html_writer, kegg)
    regress.FromDatabase()
    estimators['NIST regression'] = regress
    
    if False:
        pmaps = {}
        pmaps['succinyl-CoA'] = regress.cid2pmap(91)
        pmaps['acetoacetate'] = regress.cid2pmap(164)
        pmaps['acetoacetyl-CoA'] = regress.cid2pmap(332)
        pmaps['succinate'] = regress.cid2pmap(42)
        
        for name, pmap in pmaps.iteritems():
            print name + "\n", pmap,
            print "dG'0 = ", pmap.Transform(pH=7, pMg=0, I=0.1, T=300)
            print '-'*50
    
        sys.exit(0)
        
    gc = GroupContribution(db, html_writer, kegg)
    gc.override_gc_with_measurements = True
    gc.init()
    estimators['Milo Group Contribution'] = gc
    
    for key, thermodynamics in estimators.iteritems():
        logging.info('Writing the NIST report for %s' % key)
        html_writer.write('<p><b>%s </b>\n' % key)
        html_writer.insert_toggle(key)
        html_writer.start_div(key)
        num_estimations, rmse = nist.verify_results(thermodynamics, T_range=(298, 314))
        html_writer.end_div()
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))

if __name__ == '__main__':
    main()