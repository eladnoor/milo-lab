#!/usr/bin/python

"""
    This script analyses the predictions of all the different estimation methods:
    Alberty, Hatzimanikatis, and the Milo lab Group Contribution method
"""

import logging
from toolbox.html_writer import HtmlWriter
from pygibbs.nist import Nist
from pygibbs.kegg import Kegg
from pygibbs.feist_ecoli import Feist
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.kegg_errors import KeggParseException
from toolbox.molecule import OpenBabelError
from pygibbs.thermodynamic_estimators import LoadAllEstimators


################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    estimators = LoadAllEstimators()
    html_writer = HtmlWriter("../res/nist/report.html")
    nist = Nist()
    nist.T_range = (273.15 + 24, 273.15 + 40)
    #nist.override_I = 0.25
    #nist.override_pMg = 14.0
    #nist.override_T = 298.15
    
    html_writer.write('<p>\n')
    html_writer.write("Total number of reaction in NIST: %d</br>\n" % len(nist.data))
    html_writer.write("Total number of reaction in range %.1fK < T < %.1fK: %d</br>\n" % \
                      (nist.T_range[0], nist.T_range[1], len(nist.SelectRowsFromNist())))
    html_writer.write('</p>\n')

    reactions = {}
    reactions['KEGG'] = []
    for reaction in Kegg.getInstance().AllReactions():
        try:
            reaction.Balance(balance_water=True, exception_if_unknown=True)
            reactions['KEGG'].append(reaction)
        except (KeggReactionNotBalancedException, KeggParseException, OpenBabelError):
            pass
        
    reactions['FEIST'] = Feist.FromFiles().reactions
    reactions['NIST'] = nist.GetUniqueReactionSet()
    
    pairs = [] #[('hatzi_gc', 'PGC'), ('PGC', 'PRC')] # + [('alberty', 'PRC')]
    for t1, t2 in pairs:
        logging.info('Writing the NIST report for %s vs. %s' % 
                     (estimators[t1].name, estimators[t2].name))
        html_writer.write('<p><b>%s vs. %s</b> ' % 
                     (estimators[t1].name, estimators[t2].name))
        html_writer.insert_toggle(start_here=True)
        nist.two_way_comparison(html_writer=html_writer, 
                                thermo1=estimators[t1],
                                thermo2=estimators[t2],
                                name='%s_vs_%s' % (t1, t2))
        html_writer.div_end()
        html_writer.write('</p>')
    
    if False:
        estimators['alberty'].CompareOverKegg(html_writer, 
                                              other=estimators['PRC'],
                                              fig_name='kegg_compare_alberty_vs_nist')
    
    rowdicts = []
    rowdict = {'Method': 'Total'}
    for db_name, reaction_list in reactions.iteritems():
        rowdict[db_name + ' coverage'] = len(reaction_list)
    rowdicts.append(rowdict)
    
    for name in ['merged', 'alberty', 'PGC', 'PRC', 'hatzi_gc']:
        thermo = estimators[name]
        logging.info('Writing the NIST report for %s' % thermo.name)
        html_writer.write('<p><b>%s</b> ' % thermo.name)
        html_writer.insert_toggle(start_here=True)
        num_estimations, rmse = nist.verify_results(html_writer=html_writer, 
                                                    thermodynamics=thermo,
                                                    name=name)
        html_writer.div_end()
        html_writer.write('N = %d, RMSE = %.1f</p>\n' % (num_estimations, rmse))
        logging.info('N = %d, RMSE = %.1f' % (num_estimations, rmse))
        
        rowdict = {'Method':thermo.name,
            'RMSE (kJ/mol)':"%.1f (N=%d)" % (rmse, num_estimations)}
        for db_name, reaction_list in reactions.iteritems():
            n_covered = thermo.CalculateCoverage(reaction_list)
            percent = n_covered * 100.0 / len(reaction_list)
            rowdict[db_name + " coverage"] = "%.1f%% (%d)" % (percent, n_covered)
            logging.info(db_name + " coverage = %.1f%%" % percent)
        rowdicts.append(rowdict)
    
    headers = ['Method', 'RMSE (kJ/mol)'] + \
        [db_name + ' coverage' for db_name in reactions.keys()]
    html_writer.write_table(rowdicts, headers=headers)

if __name__ == '__main__':
    main()
