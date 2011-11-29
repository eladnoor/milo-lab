#!/usr/bin/python

import copy
import logging
import matplotlib
import os
from pylab import meshgrid
import re
import sys
import numpy as np

from copy import deepcopy
from optparse import OptionParser
from pygibbs.feasibility import pC_to_range, find_mtdf, find_pCr
from pygibbs.feasibility import LinProgNoSolutionException, find_ratio
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.kegg import Kegg
from pygibbs.pathway import PathwayData
from pygibbs.thermodynamic_constants import transform
from pygibbs.thermodynamic_constants import default_T, default_pH
from pygibbs.thermodynamic_constants import default_I, R, F
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from toolbox.util import _mkdir
import scipy.io
import matplotlib.pyplot as plt
from pygibbs.pathway_modelling import KeggPathway,\
    UnsolvableConvexProblemException, DeltaGNormalization
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.compound_abundance import CompoundAbundance
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy

class ThermodynamicAnalysis(object):
    def __init__(self, db, html_writer, thermodynamics):
        self.db = db
        self.html_writer = html_writer
        self.thermo = thermodynamics
        self.kegg = Kegg.getInstance()

        # set the standard redox potential to 320mV and concentrations to 1M 
        # the formation energy will be used only for the dG in the tables
        # but will later be overridden by the value or 'redox' which is
        # determined by the Y-axis in the contour plot. 
        default_E_prime = -0.32 # the E' of NAD(P) at pH 7
        self.thermo.AddPseudoisomer(28, nH=0, z=0, nMg=0, dG0=0) # oxidized electron carrier
        self.thermo.AddPseudoisomer(30, nH=0, z=0, nMg=0, 
                                    dG0=-default_E_prime * F) # reduced electron carrier

    def analyze_pathway(self, filename,
                        insert_toggles=True,
                        write_measured_concentrations=False):
        self.html_writer.write("<h1>Thermodynamic Pathway Analysis</h1>\n")
        entry2fields_map = ParsedKeggFile.FromKeggFile(filename)
        
        for key in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[key]
            p_data = PathwayData.FromFieldMap(field_map)
            
            if p_data.skip:
                logging.info("Skipping pathway: %s", key)
                continue
            try:
                self.html_writer.write('<b>%s - %s</b>\n' % (p_data.name,
                                                               p_data.analysis_type))
                if insert_toggles:
                    self.html_writer.insert_toggle(key)
                    self.html_writer.div_start(key)
            except KeyError:
                raise Exception("Both the 'NAME' and 'TYPE' fields must be defined for each pathway")

            logging.info("analyzing pathway: " + key)

            function_dict = {'PROFILE':self.analyze_profile,
                             'PCR':self.analyze_pCr,
                             'MTDF':self.analyze_mtdf,
                             'MTDF2D':self.analyze_mtdf_2d,
                             'REDOX':self.analyze_redox3,
                             'PROTONATION':self.analyze_protonation,
                             'STANDARD':self.analyze_standard_conditions}

            if p_data.analysis_type in function_dict:
                msg = function_dict[p_data.analysis_type](key, p_data)     
            else:
                raise Exception("Unknown analysis type: " + p_data.analysis_type)
            if insert_toggles:
                self.html_writer.div_end()
            if msg is not None:
                self.html_writer.write(msg)
            self.html_writer.write('</br>\n')
        
        if write_measured_concentrations:    
            self.html_writer.write('<b>Measured concentration table:</b></br>\n')
            if insert_toggles:
                div_id = self.html_writer.insert_toggle()
                self.html_writer.div_start(div_id)
            self.db.Query2HTML(self.html_writer,
                               "SELECT cid, media, 1000*concentration from compound_abundance ORDER BY cid, media",
                               column_names=["cid", "media", "concentration [mM]"])
            if insert_toggles:
                self.html_writer.div_end()

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if (len(tokens) == 0):
            return default_value
        if (len(tokens) > 1):
            raise Exception("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    def get_bounds(self, module_name, pathway_data):
        cid2bounds = {1: (1, 1)} # the default for H2O is 1
        field_map = pathway_data.field_map
        if "BOUND" in field_map:
            for line in field_map["BOUND"].strip().split('\t'):
                tokens = line.split(None)
                cid = int(tokens[0][1:])
                try:
                    b_lower = float(tokens[1])
                except ValueError:
                    b_lower = None    

                if len(tokens) == 2:
                    b_upper = b_lower
                elif len(tokens) == 3:
                    try:
                        b_upper = float(tokens[2])
                    except ValueError:
                        b_upper = None    
                else:
                    raise ValueError("Parsing error in BOUND definition for %s: %s" % \
                                     (module_name, line))
                cid2bounds[cid] = (b_lower, b_upper)
        else:
            abundance = CompoundAbundance.LoadConcentrationsFromSauer()
            for cid, bounds in abundance.GetAllBounds('glucose'):
                cid2bounds[cid] = bounds
        return cid2bounds
                
    def write_bounds_to_html(self, cid2bounds, c_range):
        self.html_writer.write("Concentration bounds:</br>\n")

        l = ["<b>Default</b>: %g M < concentration < %g M</br>\n" % (c_range)]
        for cid in sorted(cid2bounds.keys()):
            (b_lower, b_upper) = cid2bounds[cid]
            s = ""
            if b_lower == None:
                s += "-inf"
            else:
                s += "%g M" % b_lower
            s += ' < [<a href="%s">%s</a>] < ' % (self.kegg.cid2link(cid), self.kegg.cid2name(cid))
            if b_lower == None:
                s += "inf"
            else:
                s += "%g M" % b_upper
            l.append(s)
        self.html_writer.write_ul(l)

    def get_reactions(self, module_name, pathway_data):
        """
            read the list of reactions from the command file
        """
        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        cid_mapping = pathway_data.cid_mapping
        field_map = pathway_data.field_map
        
        mid = pathway_data.kegg_module_id
        if mid is not None:
            S, rids, fluxes, cids = self.kegg.get_module(pathway_data.kegg_module_id)
            for i, cid in enumerate(list(cids)):
                if cid in cid_mapping:
                    new_cid, coeff = cid_mapping[cid]
                    cids[i] = new_cid
                    S[:, i] *= coeff
            self.html_writer.write('<h3>Module <a href=http://www.genome.jp/dbget-bin/www_bget?M%05d>M%05d</a></h3>\n' % (mid, mid))       
        else:
            S, rids, fluxes, cids = self.kegg.parse_explicit_module(field_map, cid_mapping) 
        
        return S, rids, fluxes, cids

    def write_reactions_to_html(self, S, rids, fluxes, cids, show_cids=True):
        self.thermo.pH = 7
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        
        self.html_writer.write("<li>Reactions:</br><ul>\n")
        
        for r in range(S.shape[0]):
            self.html_writer.write('<li><a href=' + self.kegg.rid2link(rids[r]) + 
                                   '>%s ' % self.kegg.rid2string(rids[r]) + '</a>')
            self.html_writer.write("[x%g, &#x394;<sub>r</sub>G'<sup>0</sup> = %.1f] : " % (fluxes[r], dG0_r[r, 0]))
            self.html_writer.write(self.kegg.vector_to_hypertext(S[r, :].flat, cids, show_cids=show_cids))
            self.html_writer.write('</li>\n')
        
        v_total = np.dot(np.matrix(fluxes), S).flat
        dG0_total = np.dot(np.matrix(fluxes), dG0_r)[0,0]
        self.html_writer.write('<li><b>Total </b>')
        self.html_writer.write('[&#x394;<sub>r</sub>G<sup>0</sup> = %.1f kJ/mol] : \n' % dG0_total)
        self.html_writer.write(self.kegg.vector_to_hypertext(v_total, cids, show_cids=show_cids))
        self.html_writer.write("</li></ul></li>\n")
        
    def write_metabolic_graph(self, name, S, rids, cids):
        """
            draw a graph representation of the pathway
        """        
        Gdot = self.kegg.draw_pathway(S, rids, cids)
        Gdot.set_size("6, 8")
        #Gdot.set_dpi("90")
        self.html_writer.embed_dot(Gdot, name, width=600, height=800)
        self.html_writer.write('</br>')

    def write_formation_energies_to_html(self, cids):
        dG0_f_prime = self.thermo.GetTransformedFormationEnergies(cids)
        dict_list = []
        for i, cid in enumerate(cids):
            d = {}
            d['CID'] = '<a href="%s">C%05d</a>' % (self.kegg.cid2link(cid), cid)
            d['name'] = self.kegg.cid2name(cid)
            if np.isnan(dG0_f_prime[i, 0]):
                d["&#x394;<sub>f</sub>G'<sup>0</sup>"] = 'N/A'
            else:
                d["&#x394;<sub>f</sub>G'<sup>0</sup>"] = '%.1f' % dG0_f_prime[i, 0]
            dict_list.append(d)
        self.html_writer.write_table(dict_list, 
                headers=['CID', 'name', "&#x394;<sub>f</sub>G'<sup>0</sup>"])

    def get_conditions(self, pathway_data):
        self.thermo.pH = pathway_data.pH or self.thermo.pH
        self.thermo.I = pathway_data.I or self.thermo.I
        self.thermo.T = pathway_data.T or self.thermo.T
        self.thermo.pMg = pathway_data.pMg or self.thermo.pMg
        self.thermo.c_range = pathway_data.c_range or tuple(self.thermo.c_range)
        self.thermo.c_mid = pathway_data.c_mid or self.thermo.c_mid
        
        self.html_writer.write('Parameters:</br>\n')
        condition_list = ['pH = %g' % self.thermo.pH,
                          'Ionic strength = %g M' % self.thermo.I,
                          'pMg = %g' % self.thermo.pMg,
                          'Temperature = %g K' % self.thermo.T,
                          'Concentration range = %g - %g M' % self.thermo.c_range,
                          'Default concentration = %g M' % self.thermo.c_mid]
        self.html_writer.write_ul(condition_list)

    def analyze_profile(self, key, pathway_data):
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # read the list of conditions from the command file
        for condition in pathway_data.conditions:
            media, pH, I = condition.media, condition.pH, condition.I
            T, c0 = condition.T, condition.c0
            self.html_writer.write('<li>Conditions: media = %s, pH = %g, I = %g M, T = %g K, c0 = %g</li>\n' % (media, pH, I, T, c0))
        self.html_writer.write('</ol></li>\n')
        
        # prepare the legend for the profile graph
        legend = []
        dG_profiles = {}
        params_list = []
        for condition in pathway_data.conditions:
            for method in pathway_data.dG_methods:
                media, pH, I = condition.media, condition.pH, condition.I
                T, c0 = condition.T, condition.c0
                plot_key = method + ' dG (media=%s,pH=%g,I=%g,T=%g,c0=%g)' % (str(media), pH, I, T, c0)
                legend.append(plot_key)
                dG_profiles[plot_key] = []
                params_list.append((method, media, pH, I, T, c0, plot_key))

        (S, rids, fluxes, cids) = self.get_reactions(key, pathway_data)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>')
        self.write_metabolic_graph(key, S, rids, cids)
        
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = {}
        dG0_r = {}
        dG_f = {}
        dG_r = {}
        abundance = CompoundAbundance.LoadConcentrationsFromBennett()
        for (method, media, pH, I, T, c0, plot_key) in params_list:
            dG0_f[plot_key] = np.zeros((Nc, 1))
            dG_f[plot_key] = np.zeros((Nc, 1))
            for c in range(Nc):
                if method == "MILO":
                    dG0_f[plot_key][c] = self.thermo.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                elif method == "HATZI":
                    dG0_f[plot_key][c] = self.thermo.hatzi.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                else:
                    raise Exception("Unknown dG evaluation method: " + method)
                # add the effect of the concentration on the dG_f (from dG0_f to dG_f)
                dG_f[plot_key][c] = dG0_f[plot_key][c] + R * T * np.log(abundance.GetConcentration(cids[c], c0, media))
            dG0_r[plot_key] = np.dot(S, dG0_f[plot_key])
            dG_r[plot_key] = np.dot(S, dG_f[plot_key])
        
        # plot the profile graph
        profile_fig = plt.figure()
        profile_fig.hold(True)
        data = np.zeros((Nr + 1, len(legend)))
        for i, label in enumerate(legend):
            for r in range(1, Nr + 1):
                data[r, i] = sum(dG_r[label][:r, 0])
        plt.plot(data, fig=profile_fig)
        plt.legend(legend, loc="lower left")

        for i, rid in enumerate(rids):
            plt.text(i + 0.5, np.mean(data[i:(i + 2), 0]), rid,
                     fontsize=6, horizontalalignment='center',
                     backgroundcolor='white', fig=profile_fig)
        
        plt.xlabel("Reaction no.")
        plt.ylabel("dG [kJ/mol]")
        self.html_writer.embed_matplotlib_figure(profile_fig, name=key+'_prfl')
    
    def analyze_pCr(self, key, pathway_data):
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc)
        c_mid = 1e-3
        pH, I, T = default_pH, default_I, default_T
        concentration_bounds = copy.deepcopy(self.kegg.cid2bounds)
        if len(pathway_data.conditions) > 1:
            raise Exception('More than 1 condition listed for pCr analysis')
        
        if pathway_data.conditions:
            c = pathway_data.conditions[0]
            pH, I, T = c.pH, c.I, c.T
            self.html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            
        if pathway_data.c_mid:
            c_mid = pathway_data.c_mid
            
        self.html_writer.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>\n')
        self.write_metabolic_graph(key, S, rids, cids)

        field_map = pathway_data.field_map        
        physiological_pC = field_map.GetFloatField('PHYSIO', default_value=4)
        Nr, Nc = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        pC = np.arange(0, 20, 0.1)
        B_vec = np.zeros(len(pC))
        #label_vec = [""] * len(pC)
        #limiting_reactions = set()
        for i in xrange(len(pC)):
            c_range = pC_to_range(pC[i], c_mid=c_mid)
            unused_dG_f, unused_concentrations, B = find_mtdf(S, dG0_f, c_range=c_range, bounds=bounds)
            B_vec[i] = B
            #curr_limiting_reactions = set(find(abs(dG_r - B) < 1e-9)).difference(limiting_reactions)
            #label_vec[i] = ", ".join(["%d" % rids[r] for r in curr_limiting_reactions]) # all RIDs of reactions that have dG_r = B
            #limiting_reactions |= curr_limiting_reactions

        try:
            unused_dG_f, unused_concentrations, pCr = find_pCr(S, dG0_f, c_mid=c_mid, bounds=bounds)
        except LinProgNoSolutionException:
            pCr = None
            
        try:
            c_range = pC_to_range(physiological_pC, c_mid=c_mid)
            unused_dG_f, unused_concentrations, B_physiological = find_mtdf(S, dG0_f, c_range=c_range, bounds=bounds) 
        except LinProgNoSolutionException:
            B_physiological = None
        
        slack_fig = plt.figure()
        plt.plot(pC, B_vec, 'b', fig=slack_fig)
        #for i in xrange(len(pC)):
        #    plt.text(pC[i], B_vec[i], label_vec[i], fontsize=6, 
        #             horizontalalignment='left', backgroundcolor='white')
        plt.xlabel('pC')
        plt.ylabel('slack [kJ/mol]')
        (ymin, _) = plt.ylim()
        (xmin, _) = plt.xlim()
#        broken_barh([(xmin, pCr), (pCr, xmax)], (ymin, 0), facecolors=('yellow', 'green'), alpha=0.3)
        plt.axhspan(ymin, 0, facecolor='b', alpha=0.15)
        title = 'C_mid = %g' % c_mid
        if (pCr != None and pCr < pC.max()):
            title += ', pCr = %.1f' % pCr
            plt.plot([pCr, pCr], [ymin, 0], 'k--')
            plt.text(pCr, 0, 'pCr = %.1f' % pCr, fontsize=8)
            if (pCr < physiological_pC):
                plt.axvspan(pCr, physiological_pC, facecolor='g', alpha=0.3)
        if (B_physiological != None and physiological_pC < pC.max()):
            title += ', slack = %.1f [kJ/mol]' % B_physiological
            plt.plot([xmin, physiological_pC], [B_physiological, B_physiological], 'k--')
            plt.text(physiological_pC, B_physiological, 'B=%.1f' % B_physiological, fontsize=8)
        
        plt.title(title)
        plt.ylim(ymin=ymin)
        self.html_writer.embed_matplotlib_figure(slack_fig, width=800, height=600)

        # write a table of the compounds and their dG0_f
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f' [kJ/mol]"))
        for c in range(Nc):
            compound = self.kegg.cid2compound(cids[c])
            cid_str = '<a href="%s">C%05d</a>' % (compound.get_link(), compound.cid)
            self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%.1f</td>\n' % (cid_str, compound.name, dG0_f[c, 0]))
        self.html_writer.write('</table><br>\n')
        
        # write a table of the reactions and their dG0_r
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG RID", "Reaction", "flux"))
        for r in range(Nr):
            rid_str = '<a href="http://www.genome.jp/dbget-bin/www_bget?rn:R%05d">R%05d</a>' % (rids[r], rids[r])
            spr = {}
            for c in np.nonzero(S[r, :])[0]:
                spr[cids[c]] = S[r, c]
            reaction_str = self.kegg.sparse_to_hypertext(spr)
            self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%g</td>\n' % (rid_str, reaction_str, fluxes[r]))
        self.html_writer.write('</table><br>\n')

    def analyze_mtdf(self, key, pathway_data):
        self.get_conditions(pathway_data)
        cid2bounds = self.get_bounds(key, pathway_data)
        #self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        
        keggpath = KeggPathway(S, rids, fluxes, cids, None, dG0_r,
                               cid2bounds=cid2bounds, c_range=self.thermo.c_range)
        try:
            min_total_dG_prime, max_total_dG_prime = keggpath.GetTotalReactionEnergy()
            
            _, concentrations, mtdf = keggpath.FindMtdf(
                                    normalization=DeltaGNormalization.SIGN_FLUX)
        except UnsolvableConvexProblemException as e:
            self.html_writer.write("<b>WARNING: cannot calculate MTDF "
                                   "because %s:</b></br>\n" %
                                   str(e))
            problem_str = str(e.problem).replace('\n', '</br>\n')
            self.html_writer.write("%s" % problem_str)
            return
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('MTDF = %.1f [kJ/mol]' % mtdf, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, name=key+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        concentration_fig = keggpath.PlotConcentrations(concentrations)
        plt.title('MTDF = %.1f [kJ/mol]' % mtdf, figure=concentration_fig)
        self.html_writer.embed_matplotlib_figure(concentration_fig, name=key+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)

        self.write_metabolic_graph(key+"_grph", S, rids, cids)
        self.write_formation_energies_to_html(cids)
        dG_r_prime = keggpath.CalculateReactionEnergiesUsingConcentrations(concentrations)
        return "MTDF = %.1f kJ/mol, Total &#x394;<sub>r</sub>G' = %.1f [min = %.1f, max = %.1f] kJ/mol" % \
            (mtdf, float(np.dot(dG_r_prime.T, fluxes)), 
             min_total_dG_prime, max_total_dG_prime)

    def analyze_protonation(self, key, pathway_data):
        field_map = pathway_data.field_map
        pH_list = pathway_data.pH_values
        I = pathway_data.I        
        T = pathway_data.T
        
        cid = field_map.GetStringField("COMPOUND")
        cid = int(cid[1:])
        
        pmatrix = self.thermo.cid2PseudoisomerMap(cid).ToMatrix()
        data = np.zeros((len(self.cid), len(pH_list)))
        for j in range(len(pH_list)):
            pH = pH_list[j]
            nMg = 0
            pMg = 14
            dG0_array = np.matrix([-transform(dG0, nH, z, nMg, pH, pMg, I, T) / (R * T) \
                                   for (nH, z, dG0) in self.cid])
            dG0_array = dG0_array - max(dG0_array)
            p_array = np.exp(dG0_array)
            p_array = p_array / sum(p_array)
            data[:, j] = p_array    
        
        protonation_fig = plt.figure()
        plt.plot(pH_list, data.T)
        prop = matplotlib.font_manager.FontProperties(size=10)
        name = self.kegg.cid2name(cid)
        plt.legend(['%s [%d]' % (name, z) for (nH, z, dG0) in pmatrix], prop=prop)
        plt.xlabel("pH")
        plt.ylabel("Pseudoisomer proportion")
        self.html_writer.embed_matplotlib_figure(protonation_fig, width=800, height=600)
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for (nH, z, dG0) in pmatrix:
            self.html_writer.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        self.html_writer.write('</table>')

    def analyze_redox(self, key, pathway_data):
        self.thermo.I = pathway_data.I or self.thermo.I
        self.thermo.T = pathway_data.T or self.thermo.T 
        pH_list = pathway_data.pH_values
        redox_list = pathway_data.redox_values or np.arange(-3.0, 3.01, 0.5)
        c_mid = pathway_data.c_mid or thermo.c_mid

        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)

        pCr_mat = np.zeros((len(pH_list), len(redox_list)))
        for i, pH in enumerate(pH_list):
            self.thermo.pH = pH
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, redox in enumerate(redox_list):
                cid2bounds = deepcopy(self.kegg.cid2bounds)
                r = 10**(redox/2.0)
                cid2bounds[3] = (c_mid * r, c_mid * r) # NAD+
                cid2bounds[4] = (c_mid / r, c_mid / r) # NADH
                cid2bounds[6] = (c_mid * r, c_mid * r) # NADPH
                cid2bounds[5] = (c_mid / r, c_mid / r) # NADP+
                bounds = [cid2bounds.get(cid, (None, None)) for cid in cids]
                try:
                    unused_dG_f, unused_concentrations, pCr = find_pCr(S, dG0_f, c_mid=c_mid, bounds=bounds)
                except LinProgNoSolutionException:
                    pCr = -1
                pCr_mat[i, j] = pCr
                
        contour_fig = plt.figure()
        pH_meshlist, r_meshlist = meshgrid(pH_list, redox_list)
        CS = plt.contour(pH_meshlist.T, r_meshlist.T, pCr_mat)       
        plt.clabel(CS, inline=1, fontsize=10)
        plt.xlabel("pH")
        plt.ylabel("$\\log{\\frac{[NAD(P)^+]}{[NAD(P)H]}}$")
        plt.title("pCr as a function of pH and Redox state")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=800, height=600)
        
    def analyze_redox2(self, key, pathway_data):
        self.thermo.I = pathway_data.I or self.thermo.I
        self.thermo.T = pathway_data.T or self.thermo.T 
        pH_list = pathway_data.pH_values
        c_range = pathway_data.c_range or tuple(self.thermo.c_range)

        field_map = pathway_data.field_map        
        co2_list = field_map.GetVFloatField("LOG_CO2", np.arange(-5.0, -1.99, 0.25))
        
        self.html_writer.write('Parameters:</br>\n')
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>ionic strength = %g M</li>\n' % self.thermo.I)
        self.html_writer.write('<li>temperature = %g K</li>\n' % self.thermo.T)
        self.html_writer.write('<li>concentration range = %g - %g M</li>\n' % \
                               (c_range[0], c_range[1]))
        self.html_writer.write('</ul>\n')
        
        self.html_writer.insert_toggle(key)
        self.html_writer.div_start(key)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        self.html_writer.div_end()
        
        cid2bounds = {}
        cid2bounds[1] = (1.0, 1.0) # H2O
        cid2bounds[2] = (1e-3, 1e-3) # ATP
        cid2bounds[8] = (1e-4, 1e-4) # ADP
        cid2bounds[20] = (None, None) # AMP
        cid2bounds[9] = (1e-3, 1e-3) # Pi
        cid2bounds[13] = (None, None) # PPi

        cid2bounds[288] = (None, None) # Co2(tot)

        cid2bounds[3] = (None, None) # NAD+
        cid2bounds[4] = (None, None) # NADH
        
        cid2bounds[5] = (None, None) # NADPH
        cid2bounds[6] = (None, None) # NADP+
        cid2bounds[139] = (None, None)  # Ferrodoxin(ox)
        cid2bounds[138] = (None, None)  # Ferrodoxin(red)
        cid2bounds[399] = (None, None)  # Ubiquinone-10(ox)
        cid2bounds[390] = (None, None)  # Ubiquinone-10(red)
        cid2bounds[828] = (None, None)  # Menaquinone(ox)
        cid2bounds[5819] = (None, None) # Menaquinone(red)
        
        ratio_mat = np.zeros((len(pH_list), len(co2_list)))
        feasability_mat = np.zeros((len(pH_list), len(co2_list)))
        for i, pH in enumerate(pH_list):
            self.thermo.pH = pH
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, log_co2 in enumerate(co2_list):
                c_co2 = 10**(log_co2)
                cid2bounds[11] = (c_co2, c_co2) # Co2(aq)
                try:
                    _, _, log_ratio = find_ratio(S, rids, fluxes, cids, dG0_f, cid_up=4, cid_down=3, 
                        c_range=c_range, T=self.thermo.T, cid2bounds=cid2bounds)
                    ratio_mat[i, j] = log_ratio
                    feasability_mat[i, j] = 1
                except LinProgNoSolutionException:
                    ratio_mat[i, j] = 10 # 10 is a very high value for the ratio
                    feasability_mat[i, j] = 0
                
        contour_fig = plt.figure()
        contour_fig.hold(True)
        pH_meshlist, atp_meshlist = meshgrid(pH_list, co2_list)
        #plt.contourf(pH_meshlist.T, atp_meshlist.T, feasability_mat, colors=('red','white'))
        #CS = plt.contour(pH_meshlist.T, atp_meshlist.T, ratio_mat)
        plt.pcolor(pH_meshlist.T, atp_meshlist.T, ratio_mat)
        plt.colorbar()
        #plt.clabel(CS, inline=1, fontsize=10)
        plt.xlabel("pH")
        plt.ylabel("$\\log{[CO_2]}$")
        plt.title("minimal $\\log{\\frac{[NADH]}{[NAD+]}}$ required for feasibility")
        self.html_writer.embed_matplotlib_figure(contour_fig, name=key+'mesh')
        
    def analyze_redox3(self, key, pathway_data):
        """
            Plot the minimal concentration of a compound (e.g. CO2) required
            for the pathway to be feasible, at different pH and redox states.
            The redox state is defined as the E' of NADP(ox) -> NADP(red).
        """
        self.get_conditions(pathway_data)
        cid2bounds = self.get_bounds(key, pathway_data)
        self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        self.html_writer.write('</br>\n')
        
        pH_list = pathway_data.pH_values or np.arange(5.0, 9.01, 0.1)
        redox_list = pathway_data.redox_values or np.arange(-0.500, -0.249999, 0.01)
        
        pH_mat = np.zeros((len(pH_list), len(redox_list)))
        redox_mat = np.zeros((len(pH_list), len(redox_list)))
        ratio_mat = np.zeros((len(pH_list), len(redox_list)))
        
        cid_to_minimize = 11 # CO2
        max_allowed_concentration = 1.0 # 1M
        cid2bounds[cid_to_minimize] = (1e-50, max_allowed_concentration)
        
        for i, pH in enumerate(pH_list):
            self.thermo.SetConditions(pH=pH)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, redox in enumerate(redox_list):
                pH_mat[i, j] = pH
                redox_mat[i, j] = redox
                
                E0 = -0.32 # for NADP+/NADPH (in V)
                r = np.exp(-2*F/(R*default_T) * (redox - E0))
                cid2bounds[6] = (1e-5, 1e-5) # NADP+
                cid2bounds[5] = (1e-5*r, 1e-5*r) # NADPH
                logging.debug("E = %g, ratio = %.1g" % (redox, r))
                
                keggpath = KeggPathway(S, rids, fluxes, cids,
                                       formation_energies=dG0_f,
                                       cid2bounds=cid2bounds,
                                       c_range=self.thermo.c_range)
                try:
                    # minimize CO2 concentration
                    _, _, min_ln_conc = keggpath.FindMinimalFeasibleConcentration(cid_to_minimize)
                    ratio_mat[i, j] = min_ln_conc / np.log(10) # change to base 10
                except UnsolvableConvexProblemException as e:
                    logging.debug(str(e))
                    ratio_mat[i, j] = np.log10(max_allowed_concentration) + \
                                      1.0 # i.e. 10 times higher than the upper bound
        
        field_map = pathway_data.field_map  
        matfile = field_map.GetStringField("MATFILE", "")
        if matfile:
            scipy.io.savemat(matfile, {"pH":pH_mat, "redox":redox_mat, 
                "co2":ratio_mat, "S":S, "rids":np.array(rids),
                "fluxes":np.array(fluxes), "cids":np.array(cids),
                "dG0_f":dG0_f}, oned_as='column')
        
        fig = plt.figure()
        fig.hold(True)
        #CS = plt.contour(pH_mat, redox_mat, ratio_mat, np.arange(-4.5, 0.01, 0.5), colors='k')
        #plt.clabel(CS, inline=1, fontsize=7, colors='black')
        collection = plt.pcolor(pH_mat, redox_mat, ratio_mat, figure=fig)
        plt.colorbar(mappable=collection)
        plt.xlim(min(pH_list), max(pH_list))
        plt.ylim(min(redox_list), max(redox_list))
        plt.xlabel("pH")
        plt.ylabel("$E^'$ (V)")
        plt.title("minimal $\\log([CO_2]/1M)$ required for feasibility")
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480)

    def analyze_mtdf_2d(self, key, pathway_data, contour=True):
        """
            Plot the MTDF of a pathway, at different pH and redox states.
            The redox state is defined as the E' of NADP(ox) -> NADP(red).
        """
        self.get_conditions(pathway_data)
        cid2bounds = self.get_bounds(key, pathway_data)

        # since the E' is one of the parameters being modulated (Y-axis)
        # there is no sense in allowing the concentrations to vary, since that
        # has an equivalent effect on the dG'.
        cid2bounds[28] = (1.0, 1.0)
        cid2bounds[30] = (1.0, 1.0)

        self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        self.html_writer.write('</br>\n')

        dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
        logging.warning('acceptor(red) = %.1f' % dG0_f[cids.index(30)])
        logging.warning('acceptor(ox) = %.1f' % dG0_f[cids.index(28)])
        logging.warning('-E*F = %.1f' % (0.32 * F))
        logging.warning('pH = %.1f' % self.thermo.pH)

        keggpath = KeggPathway(S, rids, fluxes, cids,
                               formation_energies=dG0_f,
                               cid2bounds=cid2bounds,
                               c_range=self.thermo.c_range)
        try:
            _, concentrations, mtdf = keggpath.FindMtdf()
        except UnsolvableConvexProblemException as e:
            self.html_writer.write("<b>WARNING: cannot calculate MTDF "
                                   "because %s:</b></br>\n" %
                                   str(e))
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('MTDF = %.1f [kJ/mol]' % mtdf, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, name=key+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        concentration_fig = keggpath.PlotConcentrations(concentrations)
        plt.title('MTDF = %.1f [kJ/mol]' % mtdf, figure=concentration_fig)
        self.html_writer.embed_matplotlib_figure(concentration_fig, name=key+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)
        
        self.html_writer.write('</br>\n')
        
        pH_list = pathway_data.pH_values or np.arange(5.0, 9.01, 0.25)
        redox_list = pathway_data.redox_values or np.arange(-0.500, 0.2001, 0.05)
        
        pH_mat = np.zeros((len(pH_list), len(redox_list)))
        redox_mat = np.zeros((len(pH_list), len(redox_list)))
        mtdf_mat = np.zeros((len(pH_list), len(redox_list)))
        
        cid_to_minimize = 11 # CO2
        max_allowed_concentration = 1.0 # 1M
        cid2bounds[cid_to_minimize] = (1e-50, max_allowed_concentration)
        
        for i, pH in enumerate(pH_list):
            self.thermo.SetConditions(pH=pH)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, E in enumerate(redox_list):
                pH_mat[i, j] = pH
                redox_mat[i, j] = E
                
                # override the formation energy of the reduced species to match
                # the value of E.
                dG0_f[cids.index(28)] = 0
                dG0_f[cids.index(30)] = -E * F
                keggpath = KeggPathway(S, rids, fluxes, cids,
                                       formation_energies=dG0_f,
                                       cid2bounds=cid2bounds,
                                       c_range=self.thermo.c_range)
                try:
                    _, _, mtdf = keggpath.FindMtdf()
                    mtdf_mat[i, j] = mtdf
                    
                    logging.warning('E = %g, pH = %g, MTDF = %.1f' % (E, pH, mtdf))
                    
                except UnsolvableConvexProblemException as e:
                    logging.debug(str(e))
                    mtdf_mat[i, j] = np.NaN
        
        fig = plt.figure()
        fig.hold(True)
        if contour:
            CS = plt.contour(pH_mat, redox_mat, mtdf_mat, 
                             np.arange(-2.5, 20.01, 2.5), colors='k')
            plt.clabel(CS, inline=1, fontsize=7, colors='black')
        else:
            collection = plt.pcolor(pH_mat, redox_mat, mtdf_mat, figure=fig)
            plt.colorbar(mappable=collection)
        plt.xlim(min(pH_list), max(pH_list))
        plt.ylim(min(redox_list), max(redox_list))
        plt.xlabel("pH")
        plt.ylabel("$E^'$ (V)")
        plt.title("MTDF (in kJ/mol)")
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480,
                                                 name=key+"_cntr")

    def analyze_standard_conditions(self, key, pathway_data):
        self.get_conditions(pathway_data)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)


def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="merged",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          default="../data/thermodynamics/pathways.txt",
                          help="The file to read for pathways to analyze.")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          default='../res/thermo_analysis/report.html',
                          help="Where to write output to.")
    return opt_parser


if __name__ == "__main__":
    plt.rcParams['text.usetex'] = False
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 5
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams['figure.figsize'] = [6.0, 6.0]
    plt.rcParams['figure.dpi'] = 100

    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    input_filename = os.path.abspath(options.input_filename)
    output_filename = os.path.abspath(options.output_filename)
    if not os.path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename
    print 'Will write output to %s' % output_filename
    
    db_loc = options.db_filename
    print 'Reading from DB %s' % db_loc
    db = SqliteDatabase(db_loc)

    thermo = estimators[options.thermodynamics_source]
    print "Using the thermodynamic estimations of: " + thermo.name
    
    kegg = Kegg.getInstance()
    thermo.bounds = deepcopy(kegg.cid2bounds)
    
    dirname = os.path.dirname(output_filename)
    if not os.path.exists(dirname):
        print 'Making output directory %s' % dirname
        _mkdir(dirname)
    
    print 'Executing thermodynamic pathway analysis'
    html_writer = HtmlWriter(output_filename)
    thermo_analyze = ThermodynamicAnalysis(db, html_writer, thermodynamics=thermo)
    thermo_analyze.analyze_pathway(input_filename)

    
