#!/usr/bin/python

"""
    This script analyses the predictions of all the different estimation methods:
    Alberty, Hatzimanikatis, and the Milo lab Group Contribution method
"""

from pylab import * #@UnusedWildImport
from toolbox.html_writer import HtmlWriter
from pygibbs.groups import GroupContribution
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from toolbox import database
import logging
from pygibbs.nist_regression import NistRegression
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics

################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    db_public = database.SqliteDatabase('../data/public_data.sqlite')
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter("../res/nist/report.html")
    nist = Nist(db, html_writer)
    nist.Load()
    nist.T_range = (298, 314)
    #nist.override_I = 0.25
    nist.override_pMg = 10.0

    estimators = {}
    estimators['Alberty'] = PsuedoisomerTableThermodynamics.FromDatabase(
                                db_public, 'alberty_pseudoisomers')
    estimators['Hatzimanikatis Group Contribution'] = Hatzi()

    regress = NistRegression(db, html_writer, nist) 
    regress.FromDatabase()
    estimators['NIST regression'] = regress
        
    gc = GroupContribution(db, html_writer)
    gc.override_gc_with_measurements = True
    gc.init()
    estimators['Milo Group Contribution'] = gc
    
    for key, thermodynamics in estimators.iteritems():
        logging.info('Writing the NIST report for %s' % key)
        html_writer.write('<p><b>%s </b>\n' % key)
        html_writer.insert_toggle(key)
        html_writer.start_div(key)
        num_estimations, rmse = nist.verify_results(thermodynamics)
        html_writer.end_div()
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))
        logging.info('N = %d, RMSE = %.1f' % (num_estimations, rmse))

if __name__ == '__main__':
    main()
