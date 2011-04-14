#!/usr/bin/python

"""
    This script analyses the predictions of all the different estimation methods:
    Alberty, Hatzimanikatis, and the Milo lab Group Contribution method
"""

import logging
from toolbox.html_writer import HtmlWriter
from toolbox.database import SqliteDatabase
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.kegg import Kegg

################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    db_public = SqliteDatabase('../data/public_data.sqlite')
    db_gibbs = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter("../res/nist/report.html")
    nist = Nist()
    nist.T_range = (273.15 + 24, 273.15 + 40)
    #nist.override_I = 0.25
    nist.override_pMg = 14.0
    
    html_writer.write('<p>\n')
    html_writer.write("Total number of reaction in NIST: %d</br>\n" % len(nist.data))
    html_writer.write("Total number of reaction in range %.1fK < T < %.1fK: %d</br>\n" % \
                      (nist.T_range[0], nist.T_range[1], len(nist.SelectRowsFromNist())))
    html_writer.write('</p>\n')

    db_tables = {'Alberty': (db_public, 'alberty_pseudoisomers'),
                 'NIST regression': (db_gibbs, 'nist_regression_pseudoisomers'),
                 'Milo Group Contribution': (db_gibbs, 'gc_pseudoisomers')}

    estimators = {}

    for key, (db, table_name) in db_tables.iteritems():
        if db.DoesTableExist(table_name):
            estimators[key] = PsuedoisomerTableThermodynamics.FromDatabase(db, table_name)
        else:
            logging.warning('The table %s does not exist in %s' % (table_name, str(db)))
    
    estimators['Hatzimanikatis Group Contribution'] = Hatzi(use_pKa=False)
    estimators['Hatzimanikatis Group Contribution (with pKa)'] = Hatzi(use_pKa=True)
    
    estimators['Milo Group Contribution'].override_data(estimators['Alberty'])
    
    kegg_reactions = Kegg.getInstance().AllReactions()
    nist = nist.SelectRowsFromNist()
    
    for key, thermodynamics in estimators.iteritems():
        logging.info('Writing the NIST report for %s' % key)
        html_writer.write('<p><b>%s</b> ' % key)
        html_writer.insert_toggle(key)
        html_writer.start_div(key)
        num_estimations, rmse = nist.verify_results(html_writer=html_writer, 
                                                    thermodynamics=thermodynamics)
        html_writer.end_div()
        html_writer.write('Coverage: %d out of %d reaction KEGG reactions</br>\n' % 
                          (thermodynamics.CalculateCoverage(kegg_reactions), 
                           len(kegg_reactions)))
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))
        logging.info('N = %d, RMSE = %.1f' % (num_estimations, rmse))

if __name__ == '__main__':
    main()
