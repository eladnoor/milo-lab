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
from pygibbs.thermodynamic_errors import MissingReactionEnergy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import rms_flat

################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def two_way_comparison(html_writer, thermo1, thermo2, reaction_list, name=None):
    """
        Compare the estimation errors of two different evaluation methods.
    
    Write results to HTML.
    
    Args:
        thermo1: a Thermodynamics object that provides dG estimates.
        thermo2: a Thermodynamics object that provides dG estimates.
    """
    pH, pMg, I, T = (7, 14, 0.1, 298.15)
    
    total_list = []
    
    for reaction in reaction_list:
        try:
            dG0_pred1 = reaction.PredictReactionEnergy(thermo1, pH=pH, pMg=pMg, I=I, T=T)
            dG0_pred2 = reaction.PredictReactionEnergy(thermo2, pH=pH, pMg=pMg, I=I, T=T)
        except MissingReactionEnergy:
            continue
            
        total_list.append([dG0_pred1, dG0_pred2, reaction])
    
    if not total_list:
        return 0, 0
    
    # plot the profile graph
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 8
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 2
    plt.rcParams['figure.dpi'] = 100
    
    data_mat = np.array([(x[0], x[1]) for x in total_list])
    non_nan = list(np.isfinite(data_mat.sum(1)).nonzero()[0].flat)
    
    fig2 = plt.figure(figsize=(5,5))
    plt.plot(data_mat[non_nan,0], data_mat[non_nan,1], 'b.')
    rmse = rms_flat((data_mat[non_nan,0] - data_mat[non_nan,1]).flat)
    plt.text(-50, 40, r'RMSE = %.1f [kJ/mol]' % (rmse))
    plt.xlabel(r'$\Delta G_r^\circ$ from %s [kJ/mol]' % thermo1.name)
    plt.ylabel(r'$\Delta G_r^\circ$ from %s [kJ/mol]' % thermo2.name)
    plt.plot([-200, 200], [-200, 200], 'k--')
    plt.axis([-200, 200, -200, 200])
    
    html_writer.embed_matplotlib_figure(fig2, name=name+"_eval")

    table_headers = ["#", '|diff|', "dG'0 (%s)" % thermo1.name, 
                     "dG'0 (%s)" % thermo2.name,\
                     "reaction", "rid"]
    dict_list = []
    for row in total_list:
        d = {}
        if np.isnan(row[0]) or np.isnan(row[1]):
            d["|diff|"] = 0
        else:
            d["|diff|"] = abs(row[0] - row[1])
        d["dG'0 (%s)" % thermo1.name] = row[0]
        d["dG'0 (%s)" % thermo2.name] = row[1]
        d['reaction'] = row[2].to_hypertext(show_cids=True)
        if row[2].rid is not None:
            d['rid'] = '<a href="%s">R%05d</a>' % (row[2].get_link(), row[2].rid)
        else:
            d['rid'] = ''
        dict_list.append(d)
    dict_list.sort(key=lambda d:d['|diff|'], reverse=True)
    html_writer.write_table(dict_list, table_headers, decimal=1)

def main():
    html_writer = HtmlWriter("../res/nist/report.html")
    estimators = LoadAllEstimators()
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
    
    pairs = [('hatzi_gc', 'UGC')]#, ('PGC', 'PRC')] # + [('alberty', 'PRC')]
    for t1, t2 in pairs:
        logging.info('Writing the NIST report for %s vs. %s' % 
                     (estimators[t1].name, estimators[t2].name))
        html_writer.write('<p><b>%s vs. %s</b> ' % 
                     (estimators[t1].name, estimators[t2].name))
        html_writer.insert_toggle(start_here=True)
        two_way_comparison(html_writer=html_writer, 
                           thermo1=estimators[t1],
                           thermo2=estimators[t2],
                           reaction_list=reactions['FEIST'],
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
    
    for name in ['UGC', 'PGC', 'PRC', 'alberty', 'merged', 'hatzi_gc']:
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
