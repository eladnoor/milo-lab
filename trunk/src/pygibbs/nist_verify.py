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

    db_tables = {'alberty': (db_public, 'alberty_pseudoisomers', 'Alberty'),
                 'nist_regression': (db_gibbs, 'nist_regression_pseudoisomers', 'NIST regression'),
                 'milo_gc': (db_gibbs, 'gc_pseudoisomers', 'Milo Group Contribution')}

    estimators = {}

    for key, (db, table_name, thermo_name) in db_tables.iteritems():
        if db.DoesTableExist(table_name):
            estimators[key] = PsuedoisomerTableThermodynamics.FromDatabase(
                                            db, table_name, name=thermo_name)
        else:
            logging.warning('The table %s does not exist in %s' % (table_name, str(db)))
    
    estimators['hatzi_gc'] = Hatzi(use_pKa=False)
    estimators['hatzi_gc_pka'] = Hatzi(use_pKa=True)
    
    kegg_reactions = Kegg.getInstance().AllReactions()
    nist_reactions = nist.GetUniqueReactionSet()
    
    if True:
        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators['alberty'],
                                thermo2=estimators['nist_regression'],
                                name=key)

        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators['alberty'],
                                thermo2=estimators['milo_gc'],
                                name=key)

        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators['alberty'],
                                thermo2=estimators['hatzi_gc'],
                                name=key)

        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators['hatzi_gc'],
                                thermo2=estimators['milo_gc'],
                                name=key)
        
        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators['hatzi_gc'],
                                thermo2=estimators['hatzi_gc_pka'],
                                name=key)
    
    for key, thermodynamics in estimators.iteritems():
        logging.info('Writing the NIST report for %s' % thermodynamics.name)
        html_writer.write('<p><b>%s</b> ' % thermodynamics.name)
        html_writer.insert_toggle(key)
        html_writer.start_div(key)
        num_estimations, rmse = nist.verify_results(html_writer=html_writer, 
                                                    thermodynamics=thermodynamics,
                                                    name=key)
        html_writer.end_div()
        html_writer.write('Coverage: %d out of %d KEGG reactions</br>\n' % 
                          (thermodynamics.CalculateCoverage(kegg_reactions), 
                          len(kegg_reactions)))
        html_writer.write('Coverage: %d out of %d unique NIST reactions</br>\n' % 
                          (thermodynamics.CalculateCoverage(nist_reactions), 
                          len(nist_reactions)))
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))
        logging.info('N = %d, RMSE = %.1f' % (num_estimations, rmse))

if __name__ == '__main__':
    main()
