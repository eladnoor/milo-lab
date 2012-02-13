import numpy as np
import matplotlib.pyplot as plt
import logging
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData
from toolbox.html_writer import HtmlWriter
from pygibbs.pathway_modelling import KeggPathway, DeltaGNormalization
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.thermodynamic_constants import R, default_I, default_pH,\
    default_pMg, default_T
from pygibbs.kegg_reaction import Reaction

RT = R * default_T

def pareto(kegg_file, html_writer, thermo,
           pH=default_pH, I=default_I, T=default_T, pMg=default_pMg,
           plot_profile=False, section_prefix=""):
    
    entries = kegg_file.entries()
    plot_data = np.zeros((len(entries), 5)) # ODB, ODFE, min TG, max TG, sum(fluxes)
    
    html_writer.write('<h2 id="%s_tables">Individual result tables</h1>\n' % section_prefix)
    remarks = []
    good_entries = []
    for i, entry in enumerate(entries):
        field_map = kegg_file[entry]
        p_data = PathwayData.FromFieldMap(field_map)
        
        if p_data.skip:
            logging.info("Skipping pathway: %s", entry)
            remarks.append('skipping')
            continue

        #html_writer.write('<a name="%s"></a>\n' % entry)
        html_writer.write('<h3 id="%s_%s">%s</h2>\n' % (section_prefix, entry, entry))
        thermo.SetConditions(pH=pH, I=I, T=T, pMg=pMg)

        S, rids, fluxes, cids = p_data.get_explicit_reactions()
        dG0_r_prime = thermo.GetTransfromedReactionEnergies(S, cids)
        keggpath = KeggPathway(S, rids, fluxes, cids, reaction_energies=dG0_r_prime,
                               cid2bounds=thermo.bounds, c_range=thermo.c_range)

        if np.any(np.isnan(dG0_r_prime)):
            remarks.append('NaN reaction energy')
            html_writer.write('NaN reaction energy')
            keggpath.WriteProfileToHtmlTable(html_writer)
            keggpath.WriteConcentrationsToHtmlTable(html_writer)
            continue

        #keggpath.normalization = DeltaGNormalization.TIMES_FLUX
        keggpath.normalization = DeltaGNormalization.SIGN_FLUX

        _ln_conc, odb = keggpath.FindMtdf()
        odfe = 100 * np.tanh(odb / (2*R*thermo.T))

        _ln_conc, min_tg = keggpath.GetTotalReactionEnergy(odb, maximize=False) # min TG - maximal Total dG
        ln_conc, max_tg = keggpath.GetTotalReactionEnergy(odb, maximize=True) # max TG - maximal Total dG
        concentrations = np.exp(ln_conc)
        
        good_entries.append(i)
        remarks.append('okay')
        plot_data[i, :] = [odb, odfe, min_tg, max_tg, np.sum(fluxes)]

        logging.info('%20s: ODB = %.1f [kJ/mol], maxTG = %.1f [kJ/mol]' % (entry, odb, max_tg))
        html_writer.write_ul(["ODB = %.1f [kJ/mol]" % odb,
                              "ODFE = %.1f%%" % odfe,
                              "Min Total &Delta;<sub>r</sub>G' = %.1f [kJ/mol]" % min_tg,
                              "Max Total &Delta;<sub>r</sub>G' = %.1f [kJ/mol]" % max_tg])
        keggpath.WriteProfileToHtmlTable(html_writer, concentrations)
        keggpath.WriteConcentrationsToHtmlTable(html_writer, concentrations)

    html_writer.write('<h2 id="%s_summary">Summary table</h1>\n' % section_prefix)
    dict_list = [{'Name':'<a href="#%s_%s">%s</a>' % (section_prefix, entries[i], entries[i]),
                  'ODB [kJ/mol]':'%.1f' % plot_data[i, 0],
                  'ODFE':'%.1f%%' % plot_data[i, 1],
                  'Total dG\' [kJ/mol]':'%6.1f - %6.1f' % (plot_data[i, 2], plot_data[i, 3]),
                  'sum(flux)':'%g' % plot_data[i, 4],
                  'remark':remarks[i]}
                 for i in xrange(len(entries))]
    html_writer.write_table(dict_list,
        headers=['Name', 'ODB [kJ/mol]', 'ODFE', 'Total dG\' [kJ/mol]', 'sum(flux)', 'remark'])
    
    data = plot_data[good_entries, :]
    data = data[:, (3, 0)]
    labels = [entries[i] for i in good_entries] 
    return data, labels
            
def analyze(prefix, thermo):    
    kegg_file = ParsedKeggFile.FromKeggFile('../data/thermodynamics/%s.txt' % prefix)
    html_writer = HtmlWriter('../res/%s.html' % prefix)

    co2_hydration = Reaction.FromFormula("C00011 + C00001 => C00288")
    
    I, T, pMg = 0.1, 298.15, 10
    
    #pH_vec = np.arange(5, 9.001, 0.5)
    pH_vec = np.array([7])
    co2_conc_vec = np.array([1e-5])
    data_mat = []
    
    for pH in pH_vec:
        co2_hydration_dG0_prime = thermo.GetTransfromedKeggReactionEnergies([co2_hydration])[0, 0]
        for co2_conc in co2_conc_vec:
            carbonate_conc = co2_conc * np.exp(-co2_hydration_dG0_prime / (R*T))
            thermo.bounds[11] = (co2_conc, co2_conc)
            thermo.bounds[288] = (carbonate_conc, carbonate_conc)
            
            section_prefix = 'pH_%g_CO2_%g' % (pH, co2_conc*1000)
            section_title = 'pH = %g, [CO2] = %g mM' % (pH, co2_conc*1000)
            html_writer.write('<h1 id="%s_title">%s</h1>\n' %
                              (section_prefix, section_title))
            html_writer.write_ul(['<a href="#%s_tables">Individual result tables</a>' % section_prefix,
                                  '<a href="#%s_summary">Summary table</a>' % section_prefix,
                                  '<a href="#%s_figure">Summary figure</a>' % section_prefix])

            data, labels = pareto(kegg_file, html_writer, thermo,
                pH=pH, I=I, T=T, pMg=pMg,
                section_prefix=section_prefix)
            data_mat.append(data)
    
    data_mat = np.array(data_mat)
    if data_mat.shape[0] == 1:
        pareto_fig = plt.figure(figsize=(6, 6), dpi=90)
        plt.plot(data_mat[0, :, 0], data_mat[0, :, 1], '.', figure=pareto_fig)
        for i in xrange(data_mat.shape[1]):
            if data[i, 1] < 0:
                color = 'grey'
            else:
                color = 'black'
            plt.text(data_mat[0, i, 0], data_mat[0, i, 1], labels[i],
                     ha='left', va='bottom',
                     fontsize=8, color=color, figure=pareto_fig)
        plt.title(section_title, figure=pareto_fig)
    else:
        pareto_fig = plt.figure(figsize=(10, 10), dpi=90)
        for i in xrange(data_mat.shape[1]):
            plt.plot(data_mat[:, i, 0], data_mat[:, i, 1], '-', figure=pareto_fig)
            plt.text(data_mat[0, i, 0], data_mat[0, i, 1], '%g' % pH_vec[0],
                     ha='center', fontsize=6, color='black', figure=pareto_fig)
            plt.text(data_mat[-1, i, 0], data_mat[-1, i, 1], '%g' % pH_vec[-1],
                     ha='center', fontsize=6, color='black', figure=pareto_fig)
        plt.legend(labels, loc='upper right')
        plt.title('Pareto', figure=pareto_fig)
    
    plt.xlabel('Optimal Energetic Efficiency [kJ/mol]', figure=pareto_fig)
    plt.ylabel('Optimized Distributed Bottleneck [kJ/mol]', figure=pareto_fig)
    html_writer.write('<h2 id="%s_figure">Summary figure</h1>\n' % section_prefix)

    # plot the Pareto figure showing all values (including infeasible)
    html_writer.embed_matplotlib_figure(pareto_fig, name=prefix + '_0')

    # set axes to hide infeasible pathways and focus on feasible ones
    pareto_fig.axes[0].set_xlim(None, 0)
    pareto_fig.axes[0].set_ylim(0, None)
    html_writer.embed_matplotlib_figure(pareto_fig, name=prefix + '_1')
    
    html_writer.close()

if __name__ == "__main__":
    plt.rcParams['legend.fontsize'] = 6
    estimators = LoadAllEstimators()
    
    experiments = [('odb_vs_otg_reductive_nature', 'merged'),
                   ('odb_vs_otg_reductive_c1', 'merged_C1'),
                   ('odb_vs_otg_oxidative', 'merged'),
                   ('odb_vs_otg_reductive_all', 'merged')]

    for prefix, thermo_name in experiments:
        thermo = estimators[thermo_name]
        analyze(prefix, thermo)

    #analyze('../res/pathologic/GA3P => PYR/kegg_pathway.txt',
    #        '../data/thermodynamics/concentration_bounds.csv',
    #        '../res/odb_vs_otg_gap2pyr.html',
    #        estimators)