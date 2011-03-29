#!/usr/bin/python

"""
    This script analyses the predictions of all the different estimation methods:
    Alberty, Hatzimanikatis, and the Milo lab Group Contribution method
"""

import logging
from toolbox.html_writer import HtmlWriter
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from pygibbs.nist_regression import NistRegression
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics

################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    db_public = SqliteDatabase('../data/public_data.sqlite')
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter("../res/nist/report.html")
    nist = Nist()
    nist.T_range = (298, 314)
    #nist.override_I = 0.25
    nist.override_pMg = 10.0

    estimators = {}
    
    estimators['Alberty'] = PsuedoisomerTableThermodynamics.FromDatabase(
                                db_public, 'alberty_pseudoisomers')
    
    estimators['Hatzimanikatis Group Contribution'] = Hatzi()
    estimators['Hatzimanikatis Group Contribution'].use_pKa = False

    estimators['Hatzimanikatis Group Contribution (with pKa)'] = Hatzi()
    estimators['Hatzimanikatis Group Contribution (with pKa)'].use_pKa = True
    
    estimators['NIST regression'] = NistRegression(db, html_writer, nist=nist) 
    estimators['NIST regression'].FromDatabase()
        
    estimators['Milo Group Contribution'] = GroupContribution(db, html_writer)
    estimators['Milo Group Contribution'].override_gc_with_measurements = True
    estimators['Milo Group Contribution'].init()
    
    for key, thermodynamics in estimators.iteritems():
        logging.info('Writing the NIST report for %s' % key)
        html_writer.write('<p><b>%s </b>\n' % key)
        html_writer.insert_toggle(key)
        html_writer.start_div(key)
        num_estimations, rmse = nist.verify_results(html_writer=html_writer, 
                                                    thermodynamics=thermodynamics)
        html_writer.end_div()
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))
        logging.info('N = %d, RMSE = %.1f' % (num_estimations, rmse))

if __name__ == '__main__':
    main()
