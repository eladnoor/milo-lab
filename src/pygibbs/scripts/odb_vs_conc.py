import numpy as np
import matplotlib.pyplot as plt
import logging
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from pygibbs.pathway_modelling import KeggPathway, DeltaGNormalization
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.thermodynamic_constants import R, symbol_dr_G_prime 
from argparse import ArgumentParser
import csv
from toolbox import color
from collections import defaultdict

def KeggFile2PathwayList(pathway_file):
    kegg_file = ParsedKeggFile.FromKeggFile(pathway_file)
    entries = kegg_file.entries()
    pathway_list = []
    for entry in entries:
        p_data = PathwayData.FromFieldMap(kegg_file[entry])
        pathway_list.append((entry, p_data))
    return pathway_list

def Pareto(pathway_list, html_writer, thermo, pH=None,
           plot_profile=False, section_prefix="", balance_water=True,
           override_bounds={}):
    """
        Return values are (data, labels)
        
        data - a matrix where the rows are the pathways (entries) and the columns are
               max Total dG and ODB
        labels - a list of the names of the pathways
    """
    
    
    html_writer.write('<h2 id="%s_tables">Individual result tables</h1>\n' % section_prefix)
    rowdicts = []
    for entry, p_data in pathway_list:
        rowdict = defaultdict(float) 
        rowdicts.append(rowdict)
        rowdict['entry'] = entry
        rowdict['remark'] = 'okay'

        if p_data.skip:
            logging.info("Skipping pathway: %s", rowdict['entry'])
            rowdict['remark'] = 'skipping'
            continue
        
        if pH is None:
            pH = p_data.pH
        thermo.SetConditions(pH=pH, I=p_data.I, T=p_data.T, pMg=p_data.pMg)
        thermo.c_range = p_data.c_range
        rowdict['pH'] = thermo.pH
        rowdict['I'] = thermo.I
        rowdict['T'] = thermo.T
        rowdict['pMg'] = thermo.pMg

        #html_writer.write('<a name="%s"></a>\n' % entry)
        html_writer.write('<h3 id="%s_%s">%s</h2>\n' % (section_prefix, rowdict['entry'], rowdict['entry']))

        S, rids, fluxes, cids = p_data.get_explicit_reactions(balance_water=balance_water)
        thermo.bounds = p_data.GetBounds().GetOldStyleBounds(cids)
        for cid, (lb, ub) in override_bounds.iteritems():
            thermo.bounds[cid] = (lb, ub)
        
        fluxes = np.matrix(fluxes)
        dG0_r_prime = thermo.GetTransfromedReactionEnergies(S, cids)
        keggpath = KeggPathway(S, rids, fluxes, cids, reaction_energies=dG0_r_prime,
                               cid2bounds=thermo.bounds, c_range=thermo.c_range)

        if np.any(np.isnan(dG0_r_prime)):
            html_writer.write('NaN reaction energy')
            keggpath.WriteProfileToHtmlTable(html_writer)
            keggpath.WriteConcentrationsToHtmlTable(html_writer)
            logging.info('%20s: ODB = NaN, maxTG = NaN' % (entry))
            rowdict['remark'] = 'NaN reaction energy'
            continue

        #keggpath.normalization = DeltaGNormalization.TIMES_FLUX
        keggpath.normalization = DeltaGNormalization.SIGN_FLUX

        _ln_conc, odb = keggpath.FindMtdf()
        odfe = 100 * np.tanh(odb / (2*R*thermo.T))

        _ln_conc, min_tg = keggpath.GetTotalReactionEnergy(odb, maximize=False) # min TG - maximal Total dG
        ln_conc, max_tg = keggpath.GetTotalReactionEnergy(odb, maximize=True) # max TG - maximal Total dG
        concentrations = np.exp(ln_conc)
        
        rowdict['ODB'] = odb
        rowdict['ODFE'] = odfe
        rowdict['min total dG'] = min_tg
        rowdict['max total dG'] = max_tg
        rowdict['sum of fluxes'] = np.sum(fluxes)

        logging.info('%20s: ODB = %.1f [kJ/mol], maxTG = %.1f [kJ/mol]' % (rowdict['entry'], odb, max_tg))
        html_writer.write_ul(["pH = %.1f, I = %.2fM, T = %.2f K" % (thermo.pH, thermo.I, thermo.T),
                              "ODB = %.1f [kJ/mol]" % odb,
                              "ODFE = %.1f%%" % odfe,
                              "Min Total %s = %.1f [kJ/mol]" % (symbol_dr_G_prime, min_tg),
                              "Max Total %s = %.1f [kJ/mol]" % (symbol_dr_G_prime, max_tg)])
        keggpath.WriteProfileToHtmlTable(html_writer, concentrations)
        keggpath.WriteConcentrationsToHtmlTable(html_writer, concentrations)

    html_writer.write('<h2 id="%s_summary">Summary table</h1>\n' % section_prefix)
    dict_list = [{'Name':'<a href="#%s_%s">%s</a>' % (section_prefix, d['entry'], d['entry']),
                  'ODB [kJ/mol]':'%.1f' % d['ODB'],
                  'ODFE':'%.1f%%' % d['ODFE'],
                  'Total dG\' [kJ/mol]':'%6.1f - %6.1f' % (d['min total dG'], d['max total dG']),
                  'sum(flux)':'%g' % d['sum of fluxes'],
                  'remark': d['remark']}
                 for d in rowdicts]
    html_writer.write_table(dict_list,
        headers=['Name', 'ODB [kJ/mol]', 'ODFE', 'Total dG\' [kJ/mol]', 'sum(flux)', 'remark'])

    return rowdicts
            
def analyze_conc_gradient(pathway_file, output_prefix, thermo, cids=[], pH=None):
    compound_names = ','.join([thermo.kegg.cid2name(cid) for cid in cids])
    pathway_list = KeggFile2PathwayList(pathway_file)
    pathway_names = [entry for (entry, _) in pathway_list]
    
    html_writer = HtmlWriter('%s.html' % output_prefix)
    
    # run once just to make sure that the pathways are all working:
    logging.info("testing all pathways with default concentrations")
    data = Pareto(pathway_list, html_writer, thermo,
        pH=pH, section_prefix="test", balance_water=True,
        override_bounds={})
    
    null_html_writer = NullHtmlWriter()
    csv_output = csv.writer(open('%s.csv' % output_prefix, 'w'))
    csv_output.writerow(['pH', '[' + compound_names + ']'] + pathway_names)

    conc_vec = 10**(-np.arange(2, 6.0001, 0.5)) # logarithmic scale between 10mM and 1nM
    override_bounds = {}
    
    odb_mat = []
    for conc in conc_vec.flat:
        for cid in cids:
            override_bounds[cid] = (conc, conc)
        logging.info("[%s] = %.1e M" % (compound_names, conc))
        data = Pareto(pathway_list, null_html_writer, thermo,
            pH=pH, section_prefix="", balance_water=True,
            override_bounds=override_bounds)
        odbs = [d['ODB'] for d in data]
        odb_mat.append(odbs)
        csv_output.writerow([data[0]['pH'], conc] + odbs)
    odb_mat = np.matrix(odb_mat) # rows are pathways and columns are concentrations

    fig = plt.figure(figsize=(6, 6), dpi=90)
    colormap = color.ColorMap(pathway_names)
    for i, name in enumerate(pathway_names):
        plt.plot(conc_vec, odb_mat[:, i], '.-', color=colormap[name], 
                 figure=fig)
    plt.title("ODB vs. [%s]" % (compound_names), figure=fig)
    plt.xscale('log')
    plt.xlabel('[%s] (in M)' % compound_names, figure=fig)
    plt.ylabel('Optimized Distributed Bottleneck [kJ/mol]', figure=fig)
    plt.legend(pathway_names)
    html_writer.write('<h2>Summary figure</h1>\n')
    html_writer.embed_matplotlib_figure(fig, name=output_prefix)
    
    html_writer.close()

def analyze_pH_gradient(pathway_file, output_prefix, thermo):
    pathway_list = KeggFile2PathwayList(pathway_file)
    pathway_names = [entry for (entry, _) in pathway_list]
    
    html_writer = HtmlWriter('%s.html' % output_prefix)
    
    # run once just to make sure that the pathways are all working:
    logging.info("testing all pathways with default pH")
    data = Pareto(pathway_list, html_writer, thermo,
        pH=None, section_prefix="test", balance_water=True,
        override_bounds={})
    
    null_html_writer = NullHtmlWriter()
    csv_output = csv.writer(open('%s.csv' % output_prefix, 'w'))
    csv_output.writerow(['pH'] + pathway_names)

    pH_vec = np.arange(4, 10, 0.5)
    odb_mat = []
    for pH in pH_vec.flat:
        logging.info("pH = %.1f" % (pH))
        data = Pareto(pathway_list, null_html_writer, thermo,
            pH=pH, section_prefix="", balance_water=True)
        odbs = [d['ODB'] for d in data]
        odb_mat.append(odbs)
        csv_output.writerow([data[0]['pH']] + odbs)
    odb_mat = np.matrix(odb_mat) # rows are pathways and columns are concentrations

    fig = plt.figure(figsize=(6, 6), dpi=90)
    colormap = color.ColorMap(pathway_names)
    for i, name in enumerate(pathway_names):
        plt.plot(pH_vec, odb_mat[:, i], '.-', color=colormap[name], 
                 figure=fig)
    plt.title("ODB vs. pH", figure=fig)
    plt.xlabel('pH', figure=fig)
    plt.ylabel('Optimized Distributed Bottleneck [kJ/mol]', figure=fig)
    plt.legend(pathway_names)
    html_writer.write('<h2>Summary figure</h1>\n')
    html_writer.embed_matplotlib_figure(fig, name=output_prefix)
    
    html_writer.close()

def MakeArgParser(estimators):
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser(description='A script for running ODB analysis for a gradient of concentrations')
    parser.add_argument('-c', '--cids', action='append', type=int, required=True,
                        help='an integer for the accumulator')
    parser.add_argument('-i', '--pathway_file', action='store', type=file, required=True, 
                        help="input pathway description file in KEGG format")
    parser.add_argument('-o', '--output_prefix', action='store', required=False,
                        default='../res/odb_analysis',
                        help='the prefix for the output files (*.html and *.csv)')
    parser.add_argument('-s', '--thermodynamics_source', action='store', required=False,
                        choices=estimators.keys(),
                        default="UGC",
                        help="The thermodynamic data to use")
    parser.add_argument('-p', '--pH', action="store", type=float,
                        default=None,
                        help="The pH (if provided, overrides the ones in the input file)")
    return parser

if __name__ == "__main__":
    estimators = LoadAllEstimators()
    parser = MakeArgParser(estimators)
    args = parser.parse_args()

    plt.rcParams['legend.fontsize'] = 6
    thermo = estimators[args.thermodynamics_source]

    if 80 in args.cids:
        analyze_pH_gradient(args.pathway_file, args.output_prefix, thermo)
    else:
        analyze_conc_gradient(args.pathway_file, args.output_prefix, thermo,
                              cids=args.cids, pH=args.pH)

