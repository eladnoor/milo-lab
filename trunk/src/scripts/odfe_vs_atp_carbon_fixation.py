import numpy as np
import matplotlib.pyplot as plt
import logging
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData
from toolbox.html_writer import HtmlWriter
from pygibbs.pathway_modelling import KeggPathway, DeltaGNormalization
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.metabolic_modelling.mtdf_optimizer import UnsolvableConvexProblemException
from pygibbs.thermodynamic_constants import R

def pareto(html_writer, thermo,
           filename='../data/thermodynamics/odfe_vs_atp_carbon_fixation.txt'):
    
    entry2fields_map = ParsedKeggFile.FromKeggFile(filename)
    
    c_range = (1e-6, 1e-2)
    cid2bounds = {
                  1:   (1,    1),    # water
                  11:  (1e-5, 1e-5), # CO2
                  288: (9e-5, 9e-5), # carbonate
                  2:   (5e-3, 5e-3), # ATP (in order to keep ATP -> ADP at -55 kJ/mol)
                  8:   (5e-4, 5e-4), # ADP (in order to keep ATP -> ADP at -55 kJ/mol)
                  9:   (5e-3, 5e-3), # Pi  (in order to keep ATP -> ADP at -55 kJ/mol)
                  20:  (1e-4, 1e-4), # AMP (in order to keep ATP -> AMP at -110 kJ/mol)
                  13:  (3e-9, 3e-9), # PPi (in order to keep ATP -> AMP at -110 kJ/mol)
                  #20:  (5e-3, 5e-3), # AMP (in order to keep ATP -> AMP at -55 kJ/mol)
                  #13:  (5e-2, 5e-2), # PPi (in order to keep ATP -> AMP at -55 kJ/mol)
                  3:   (1e-3, 1e-3), # NAD_ox
                  4:   (1e-4, 1e-2), # NAD_red (to keep ox/red between 0.1 and 10)
                  6:   (1e-3, 1e-3), # NADP_ox
                  5:   (1e-4, 1e-2), # NADP_red (to keep ox/red between 0.1 and 10)
                  399: (1e-3, 1e-3), # ubiquinone_ox
                  390: (1e-4, 1e-2), # ubiquinol_red (to keep ox/red between 0.1 and 10)
                  139: (1e-3, 1e-3), # ferredoxin_ox
                  138: (1e-4, 1e-2), # ferredoxin_red (to keep ox/red between 0.1 and 10)

                  #3:   (5e-3, 5e-3), # NAD_ox
                  #4:   (5e-5, 5e-5), # NAD_red (to keep ox/red at 100)
                  #6:   (5e-5, 5e-5), # NADP_ox
                  #5:   (5e-4, 5e-4), # NADP_red (to keep ox/red at 0.1)
                  #399: (1e-3, 1e-3), # ubiquinone_ox
                  #390: (1e-4, 1e-2), # ubiquinol_red (to keep ox/red between 0.1 and 10)
                  #139: (1e-6, 1e-6), # ferredoxin_ox
                  #138: (1e-2, 1e-2), # ferredoxin_red (to keep ox/red at 1e-4)
                  }
    
    keys = sorted(entry2fields_map.keys())
    plot_data = np.zeros((len(keys), 4))
    for i, key in enumerate(keys):
        field_map = entry2fields_map[key]
        p_data = PathwayData.FromFieldMap(field_map)
        
        if p_data.skip:
            logging.info("Skipping pathway: %s", key)
            continue

        html_writer.write('<h2>%s - %s</h2>\n' % (p_data.name, p_data.analysis_type))
    
        S, rids, fluxes, cids = p_data.get_explicit_reactions()
        dG0_r = thermo.GetTransfromedReactionEnergies(S, cids)
        keggpath = KeggPathway(S, rids, fluxes, cids, None, dG0_r,
                               cid2bounds=cid2bounds, c_range=c_range)
        try:
            _, _, mtdf = keggpath.FindMtdf(normalization=DeltaGNormalization.SIGN_FLUX)
            _, concentrations, total_dG_prime = keggpath.GetMaxReactionEnergy(mtdf)
            
        except UnsolvableConvexProblemException as e:
            html_writer.write("<b>WARNING: cannot calculate MTDF "
                                   "because %s:</b></br>\n" %
                                   str(e))
            problem_str = str(e.problem).replace('\n', '</br>\n')
            html_writer.write("%s" % problem_str)
            continue
        
        odfe = 100 * np.tanh(mtdf / (2*R*thermo.T))
        html_writer.write_ul(["MTDF = %.1f" % mtdf,
                              "ODFE = %.1f%%" % odfe,
                              "Total &#x394;<sub>r</sub>G' = %.1f kJ/mol" % total_dG_prime])
        #profile_fig = keggpath.PlotProfile(concentrations)
        #plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        #html_writer.embed_matplotlib_figure(profile_fig, name=key+"_prfl", height=400, width=400)
        keggpath.WriteProfileToHtmlTable(html_writer, concentrations)
        keggpath.WriteConcentrationsToHtmlTable(html_writer, concentrations)
            
        plot_data[i, :] = [mtdf, odfe, total_dG_prime, np.sum(fluxes)]
    
    fig = plt.figure(figsize=(6, 6), dpi=80)
    plt.plot(plot_data[:, 0], plot_data[:, 2], '.', figure=fig)
    for i, key in enumerate(keys):
        if plot_data[i, 0] < 0:
            color = 'red'
        else:
            color = 'black'
        plt.text(plot_data[i, 0], plot_data[i, 2]+3, key, ha='center', va='bottom', fontsize=10, color=color)
    plt.xlabel('Optimized Distributed Bottleneck [kJ/mol]', figure=fig)
    plt.ylabel('Optimal Energetic Cost [kJ/mol]', figure=fig)
    return fig
    
if __name__ == "__main__":
    thermo = LoadAllEstimators()['merged']
    output_filename = '../res/odfe_vs_atp_carbon_fixation.html'
    html_writer = HtmlWriter(output_filename)
    for pH in [5, 7]:
        thermo.SetConditions(pH=pH, I=0.25, T=298.15, pMg=10)

        #html_writer.write('<h1>pH = %g</h1>\n' % pH)
        html_writer.insert_toggle(start_here=True)
        fig = pareto(html_writer, thermo)
        html_writer.div_end()
        
        plt.title('pH = %g' % pH, figure=fig)
        html_writer.embed_matplotlib_figure(fig)
    
    html_writer.close()
