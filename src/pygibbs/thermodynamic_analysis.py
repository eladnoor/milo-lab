#!/usr/bin/python

import copy
import logging
import matplotlib
import os
import pylab
import re
import sys

import matplotlib.pyplot as plt
import scipy.io

from copy import deepcopy
from optparse import OptionParser
from pygibbs.feasibility import pC_to_range, find_mtdf, find_pCr
from pygibbs.feasibility import LinProgNoSolutionException, thermodynamic_pathway_analysis, find_ratio
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.kegg import Kegg
from pygibbs import kegg_utils
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.thermodynamic_constants import transform
from pygibbs.thermodynamic_constants import default_T, default_pH
from pygibbs.thermodynamic_constants import default_I, default_c0, R
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from toolbox.util import _mkdir


class ThermodynamicAnalysis(object):
    def __init__(self, db, html_writer, thermodynamics):
        self.db = db
        self.html_writer = html_writer
        self.thermo = thermodynamics
        self.kegg = Kegg.getInstance()

    def analyze_pathway(self, filename,
                        insert_toggles=True,
                        write_measured_concentrations=False):
        self.html_writer.write("<h1>Thermodynamic Pathway Analysis</h1>\n")
        entry2fields_map = ParsedKeggFile.FromKeggFile(filename)
        
        for key in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[key]
            if field_map.GetBoolField('SKIP'):
                logging.info("Skipping pathway: %s", key)
                continue
            try:
                self.html_writer.write('<p>\n')
                a_name = field_map.GetStringField('NAME')
                a_type = field_map.GetStringField('TYPE')
                self.html_writer.write('<h2>%s - %s</h2>\n' % (a_name, a_type))
                if insert_toggles:
                    self.html_writer.insert_toggle(key)
                    self.html_writer.start_div(key)
            except KeyError:
                raise Exception("Both the 'NAME' and 'TYPE' fields must be defined for each pathway")

            logging.info("analyzing pathway: " + key)

            function_dict = {'PROFILE':self.analyze_profile,
                             'SLACK':self.analyze_slack,
                             'MARGIN':self.analyze_margin,
                             'CONTOUR':self.analyze_contour,
                             'REDOX':self.analyze_redox3,
                             'PROTONATION':self.analyze_protonation,
                             'STANDARD':self.analyze_standard_conditions}

            analysis_type = field_map.GetStringField('TYPE')
            if analysis_type in function_dict:
                function_dict[analysis_type](key, field_map)     
            else:
                raise Exception("Unknown analysis type: " + analysis_type)
            if insert_toggles:
                self.html_writer.end_div()
            self.html_writer.write('</p>\n')
        
        if write_measured_concentrations:    
            self.html_writer.write('<p>\n')
            self.html_writer.write('<h2>Measured concentration table:</h2>\n')
            if insert_toggles:
                div_id = self.html_writer.insert_toggle()
                self.html_writer.start_div(div_id)
            self.db.Query2HTML(self.html_writer,
                               "SELECT cid, media, 1000*concentration from compound_abundance ORDER BY cid, media",
                               column_names=["cid", "media", "concentration [mM]"])
            if insert_toggles:
                self.html_writer.end_div()
            self.html_writer.write('</p>\n')

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if (len(tokens) == 0):
            return default_value
        if (len(tokens) > 1):
            raise Exception("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    def get_bounds(self, module_name, field_map):
        cid2bounds = {1: (1, 1)} # the default for H2O is 1
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

    def get_reactions(self, module_name, field_map):
        """
            read the list of reactions from the command file
        """
        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        cid_mapping = {}
        if "MAP_CID" in field_map:
            for line in field_map["MAP_CID"].strip().split('\t'):
                cid_before, cid_after = [int(cid[1:]) for cid in line.split(None, 1)]
                cid_mapping[cid_before] = (cid_after, 1.0)
        
        if "MODULE" in field_map:
            mid_str = field_map["MODULE"]
            if (mid_str[0] == 'M'):
                mid = int(mid_str[1:])
            else:
                mid = int(mid_str)
            S, rids, fluxes, cids = self.kegg.get_module(mid)
            for i, cid in enumerate(list(cids)):
                if cid in cid_mapping:
                    (new_cid, coeff) = cid_mapping[cid]
                    cids[i] = new_cid
                    S[:, i] *= coeff
            self.html_writer.write('<h3>Module <a href=http://www.genome.jp/dbget-bin/www_bget?M%05d>M%05d</a></h3>\n' % (mid, mid))       
        else:
            S, rids, fluxes, cids = self.kegg.parse_explicit_module(field_map, cid_mapping) 
        
        return (S, rids, fluxes, cids)

    def write_reactions_to_html(self, S, rids, fluxes, cids, show_cids=True):
        self.thermo.pH = 7
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        
        self.html_writer.write("<li>Reactions:</br><ul>\n")
        
        for r in range(S.shape[0]):
            self.html_writer.write('<li><a href=' + self.kegg.rid2link(rids[r]) + '>R%05d ' % rids[r] + '</a>')
            self.html_writer.write('[x%g, &#x394;G<sub>r</sub><sup>0</sup> = %.1f] : ' % (fluxes[r], dG0_r[r, 0]))
            self.html_writer.write(self.kegg.vector_to_hypertext(S[r, :].flat, cids, show_cids=show_cids))
            self.html_writer.write('</li>\n')
        
        v_total = pylab.dot(pylab.matrix(fluxes), S).flat
        dG0_total = pylab.dot(pylab.matrix(fluxes), dG0_r)[0,0]
        self.html_writer.write('<li><b>Total </b>')
        self.html_writer.write('[&#x394;G<sub>r</sub><sup>0</sup> = %.1f kJ/mol] : \n' % dG0_total)
        self.html_writer.write(self.kegg.vector_to_hypertext(v_total, cids, show_cids=show_cids))
        self.html_writer.write("</li></ul></li>\n")
        
        # Write the kegg-formatted pathway.
        reactions = map(self.kegg.rid2reaction, rids)
        kegg_utils.write_kegg_pathway(self.html_writer, reactions, fluxes)
        
    def write_metabolic_graph(self, name, S, rids, cids):
        """
            draw a graph representation of the pathway
        """        
        Gdot = self.kegg.draw_pathway(S, rids, cids)
        self.html_writer.embed_dot(Gdot, name, width=400, height=400)

    def get_conditions(self, field_map):
        self.thermo.pH = field_map.GetFloatField("PH", self.thermo.pH)
        self.thermo.I = field_map.GetFloatField("I", self.thermo.I)
        self.thermo.T = field_map.GetFloatField("T", self.thermo.T)
        self.thermo.pMg = field_map.GetFloatField("PMG", self.thermo.pMg)
        self.thermo.c_range = tuple(field_map.GetVFloatField("C_RANGE", self.thermo.c_range))
        self.thermo.c_mid = field_map.GetFloatField('C_MID', default_value=self.thermo.c_mid)
        
        self.html_writer.write('Parameters:</br>\n')
        condition_list = ['pH = %g' % self.thermo.pH,
                          'Ionic strength = %g M' % self.thermo.I,
                          'pMg = %g' % self.thermo.pMg,
                          'Temperature = %g K' % self.thermo.T,
                          'Concentration range = %g - %g M' % self.thermo.c_range,
                          'Default concentration = %g M' % self.thermo.c_mid]
        self.html_writer.write_ul(condition_list)

    def analyze_profile(self, key, field_map):
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # read the list of conditions from the command file
        conditions = []
        for condition in field_map["CONDITIONS"].split('\t'):
            (media, pH, I, T, c0) = (None, default_pH, default_I, default_T, default_c0)
            media = re.findall("media=([a-zA-Z_]+)", condition)[0]
            if (media == 'None'):
                media = None
            pH = ThermodynamicAnalysis.get_float_parameter(condition, "pH", default_pH)
            I = ThermodynamicAnalysis.get_float_parameter(condition, "I", default_I)
            T = ThermodynamicAnalysis.get_float_parameter(condition, "T", default_T)
            c0 = ThermodynamicAnalysis.get_float_parameter(condition, "c0", default_c0)
            conditions.append((media, pH, I, T, c0))
            self.html_writer.write('<li>Conditions: media = %s, pH = %g, I = %g M, T = %g K, c0 = %g</li>\n' % (media, pH, I, T, c0))
        self.html_writer.write('</ol></li>\n')
        
        # read the list of methods for calculating the dG
        methods = []
        if field_map.GetBoolField('MILO', default_value=True):
            methods.append('MILO')
        if field_map.GetBoolField('HATZI', default_value=False):
            methods.append('HATZI')
        
        # prepare the legend for the profile graph
        legend = []
        dG_profiles = {}
        params_list = []
        for (media, pH, I, T, c0) in conditions:
            for method in methods:
                plot_key = method + ' dG (media=%s,pH=%g,I=%g,T=%g,c0=%g)' % (str(media), pH, I, T, c0)
                legend.append(plot_key)
                dG_profiles[plot_key] = []
                params_list.append((method, media, pH, I, T, c0, plot_key))

        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>')
        self.write_metabolic_graph(key, S, rids, cids)
        
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = {}
        dG0_r = {}
        dG_f = {}
        dG_r = {}
        for (method, media, pH, I, T, c0, plot_key) in params_list:
            dG0_f[plot_key] = pylab.zeros((Nc, 1))
            dG_f[plot_key] = pylab.zeros((Nc, 1))
            for c in range(Nc):
                if method == "MILO":
                    dG0_f[plot_key][c] = self.thermo.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                elif method == "HATZI":
                    dG0_f[plot_key][c] = self.thermo.hatzi.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                else:
                    raise Exception("Unknown dG evaluation method: " + method)
                # add the effect of the concentration on the dG_f (from dG0_f to dG_f)
                dG_f[plot_key][c] = dG0_f[plot_key][c] + R * T * pylab.log(self.get_concentration(cids[c], c0, media))
            dG0_r[plot_key] = pylab.dot(S, dG0_f[plot_key])
            dG_r[plot_key] = pylab.dot(S, dG_f[plot_key])
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 2
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        profile_fig = pylab.figure()
        pylab.hold(True)
        data = pylab.zeros((Nr + 1, len(legend)))
        for i in range(len(legend)):
            for r in range(1, Nr + 1):
                data[r, i] = sum(dG_r[legend[i]][:r, 0])
        pylab.plot(data)
        legend(legend, loc="lower left")

        for i in range(len(rids)):
            pylab.text(i + 0.5, pylab.mean(data[i:(i + 2), 0]), rids[i], fontsize=6, horizontalalignment='center', backgroundcolor='white')
        
        pylab.xlabel("Reaction no.")
        pylab.ylabel("dG [kJ/mol]")
        self.html_writer.embed_matplotlib_figure(profile_fig, width=800, heigh=600)
    
    def analyze_slack(self, key, field_map):
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc) 
        c_mid = field_map.GetFloatField('C_MID', default_value=1e-3)
        (pH, I, T) = (default_pH, default_I, default_T)
        concentration_bounds = copy.deepcopy(self.kegg.cid2bounds)
        if ("CONDITIONS" in field_map):
            pH = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            I = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            T = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            self.html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                concentration_bounds[cid] = (conc, conc)
                self.html_writer.write(', [C%05d] = %g\n' % (cid, conc))
            self.html_writer.write('</li>\n')
        self.html_writer.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>\n')
        self.write_metabolic_graph(key, S, rids, cids)
        
        physiological_pC = field_map.GetFloatField('PHYSIO', default_value=4)
        Nr, Nc = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        pC = pylab.arange(0, 20, 0.1)
        B_vec = pylab.zeros(len(pC))
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

        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 5
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        slack_fig = pylab.figure()
        pylab.plot(pC, B_vec, 'b')
        #for i in xrange(len(pC)):
        #    text(pC[i], B_vec[i], label_vec[i], fontsize=6, horizontalalignment='left', backgroundcolor='white')
        pylab.xlabel('pC')
        pylab.ylabel('slack [kJ/mol]')
        (ymin, _) = pylab.ylim()
        (xmin, _) = pylab.xlim()
#        broken_barh([(xmin, pCr), (pCr, xmax)], (ymin, 0), facecolors=('yellow', 'green'), alpha=0.3)
        pylab.axhspan(ymin, 0, facecolor='b', alpha=0.15)
        title = 'C_mid = %g' % c_mid
        if (pCr != None and pCr < pC.max()):
            title += ', pCr = %.1f' % pCr
            pylab.plot([pCr, pCr], [ymin, 0], 'k--')
            pylab.text(pCr, 0, 'pCr = %.1f' % pCr, fontsize=8)
            if (pCr < physiological_pC):
                pylab.axvspan(pCr, physiological_pC, facecolor='g', alpha=0.3)
        if (B_physiological != None and physiological_pC < pC.max()):
            title += ', slack = %.1f [kJ/mol]' % B_physiological
            pylab.plot([xmin, physiological_pC], [B_physiological, B_physiological], 'k--')
            pylab.text(physiological_pC, B_physiological, 'B=%.1f' % B_physiological, fontsize=8)
        
        pylab.title(title)
        pylab.ylim(ymin=ymin)
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
            for c in pylab.find(S[r, :]):
                spr[cids[c]] = S[r, c]
            reaction_str = self.kegg.sparse_to_hypertext(spr)
            self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%g</td>\n' % (rid_str, reaction_str, fluxes[r]))
        self.html_writer.write('</table><br>\n')

    def analyze_margin(self, key, field_map):
        self.get_conditions(field_map)
        cid2bounds = self.get_bounds(key, field_map)
        self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
                
        thermodynamic_pathway_analysis(S, rids, fluxes, cids, self.thermo, self.html_writer)

    def analyze_contour(self, key, field_map):
        pH_list = field_map.GetVFloatField("PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I_list = field_map.GetVFloatField("I", [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        T = field_map.GetFloatField("T", default_T)
        most_abundant = field_map.GetBoolField("ABUNDANT", False)
        formula = field_map.GetStringField("REACTION")
        c0 = field_map.GetFloatField("C0", 1.0)
        media = field_map.GetStringField("MEDIA", "None")
        if media == "None":
            media = None
            
        sparse_reaction = self.kegg.formula_to_sparse(formula)
        dG_r = self.estimate_dG_reaction(sparse_reaction, pH_list, I_list, T, c0, media, most_abundant)
        contour_fig = pylab.figure()
        
        pH_meshlist, I_meshlist = pylab.meshgrid(pH_list, I_list)
        CS = pylab.contour(pH_meshlist.T, I_meshlist.T, dG_r)       
        pylab.clabel(CS, inline=1, fontsize=10)
        pylab.xlabel("pH")
        pylab.ylabel("Ionic Strength")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=800, height=600)
        self.html_writer.write('<br>\n' + self.kegg.sparse_to_hypertext(sparse_reaction) + '<br>\n')

    def analyze_protonation(self, key, field_map):
        pH_list = field_map.GetVFloatField("PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I = field_map.GetFloatField("I", default_I)
        T = field_map.GetFloatField("T", default_T)
        cid = field_map.GetStringField("COMPOUND")
        cid = int(cid[1:])
        pmatrix = self.thermo.cid2PseudoisomerMap(cid).ToMatrix()
        data = pylab.zeros((len(self.cid), len(pH_list)))
        for j in range(len(pH_list)):
            pH = pH_list[j]
            dG0_array = pylab.matrix([-transform(dG0, nH, z, pH, I, T) / (R * T) \
                                      for (nH, z, dG0) in self.cid])
            dG0_array = dG0_array - max(dG0_array)
            p_array = pylab.exp(dG0_array)
            p_array = p_array / sum(p_array)
            data[:, j] = p_array    
        
        protonation_fig = pylab.figure()
        pylab.plot(pH_list, data.T)
        prop = matplotlib.font_manager.FontProperties(size=10)
        name = self.kegg.cid2name(cid)
        pylab.legend(['%s [%d]' % (name, z) for (nH, z, dG0) in pmatrix], prop=prop)
        pylab.xlabel("pH")
        pylab.ylabel("Pseudoisomer proportion")
        self.html_writer.embed_matplotlib_figure(protonation_fig, width=800, height=600)
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for (nH, z, dG0) in pmatrix:
            self.html_writer.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        self.html_writer.write('</table>')

    def analyze_redox(self, key, field_map):
        self.thermo.I = field_map.GetFloatField("I", self.thermo.I)
        self.thermo.T = field_map.GetFloatField("T", self.thermo.T)
        pH_list = field_map.GetVFloatField("PH", pylab.arange(4.0, 10.01, 0.25))
        redox_list = field_map.GetVFloatField("REDOX", pylab.arange(-3.0, 3.01, 0.5))
        c_mid = field_map.GetFloatField("C_MID", thermo.c_mid)

        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)

        pCr_mat = pylab.zeros((len(pH_list), len(redox_list)))
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
                
        contour_fig = pylab.figure()
        pH_meshlist, r_meshlist = pylab.meshgrid(pH_list, redox_list)
        CS = pylab.contour(pH_meshlist.T, r_meshlist.T, pCr_mat)       
        pylab.clabel(CS, inline=1, fontsize=10)
        pylab.xlabel("pH")
        pylab.ylabel("$\\log{\\frac{[NAD(P)^+]}{[NAD(P)H]}}$")
        pylab.title("pCr as a function of pH and Redox state")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=800, height=600)
        
    def analyze_redox2(self, key, field_map):
        self.thermo.I = field_map.GetFloatField("I", self.thermo.I)
        self.thermo.T = field_map.GetFloatField("T", self.thermo.T)
        pH_list = field_map.GetVFloatField("PH", pylab.arange(5.0, 9.01, 0.25))
        co2_list = field_map.GetVFloatField("LOG_CO2", pylab.arange(-5.0, -1.99, 0.25))
        c_range = tuple(field_map.GetVFloatField("C_RANGE", self.thermo.c_range))
        
        self.html_writer.write('Parameters:</br>\n')
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>ionic strength = %g M</li>\n' % self.thermo.I)
        self.html_writer.write('<li>temperature = %g K</li>\n' % self.thermo.T)
        self.html_writer.write('<li>concentration range = %g - %g M</li>\n' % \
                               (c_range[0], c_range[1]))
        self.html_writer.write('</ul>\n')
        
        self.html_writer.insert_toggle(key)
        self.html_writer.start_div(key)
        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        self.html_writer.end_div()
        
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
        
        ratio_mat = pylab.zeros((len(pH_list), len(co2_list)))
        feasability_mat = pylab.zeros((len(pH_list), len(co2_list)))
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
                
        contour_fig = pylab.figure()
        pylab.hold(True)
        pH_meshlist, atp_meshlist = pylab.meshgrid(pH_list, co2_list)
        #pylab.contourf(pH_meshlist.T, atp_meshlist.T, feasability_mat, colors=('red','white'))
        #CS = pylab.contour(pH_meshlist.T, atp_meshlist.T, ratio_mat)
        pylab.pcolor(pH_meshlist.T, atp_meshlist.T, ratio_mat)
        pylab.colorbar()
        #pylab.clabel(CS, inline=1, fontsize=10)
        pylab.xlabel("pH")
        pylab.ylabel("$\\log{[CO_2]}$")
        pylab.title("minimal $\\log{\\frac{[NADH]}{[NAD+]}}$ required for feasibility")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=640, height=480)
        
    def analyze_redox3(self, key, field_map):
        self.html_writer.insert_toggle(key)
        self.html_writer.start_div(key)
        self.get_conditions(field_map)
        cid2bounds = self.get_bounds(key, field_map)
        self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)
        self.html_writer.end_div()
        self.html_writer.write('</br>\n')

        pH_list = field_map.GetVFloatField("PH", pylab.arange(5.0, 9.01, 0.2))
        redox_list = field_map.GetVFloatField("REDOX", pylab.arange(0.0, 3.01, 0.2))
        
        pH_mat = pylab.zeros((len(pH_list), len(redox_list)))
        redox_mat = pylab.zeros((len(pH_list), len(redox_list)))
        ratio_mat = pylab.zeros((len(pH_list), len(redox_list)))
        feasibility_mat = pylab.zeros((len(pH_list), len(redox_list)))
        for i, pH in enumerate(pH_list):
            self.thermo.pH = pH
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, redox in enumerate(redox_list):
                r = 10**(redox)
                cid2bounds[6] = (1e-5, 1e-5) # NADP+
                cid2bounds[5] = (1e-5*r, 1e-5*r) # NADPH

                pH_mat[i, j] = pH
                redox_mat[i, j] = redox
                
                try:
                    _, _, log_ratio = find_ratio(S, rids, fluxes, cids, dG0_f, cid_up=11, cid_down=1, 
                        c_range=self.thermo.c_range, T=self.thermo.T, cid2bounds=cid2bounds)
                    ratio_mat[i, j] = log_ratio
                    if log_ratio < -5:
                        feasibility_mat[i, j] = 0.0
                    else:
                        feasibility_mat[i, j] = 0.5
                except LinProgNoSolutionException:
                    ratio_mat[i, j] = cid2bounds[11][1] + 1.0 # i.e. 10 times higher than the upper bound
                    feasibility_mat[i, j] = 1.0
                
        matfile = field_map.GetStringField("MATFILE", "")
        if matfile:
            scipy.io.savemat(matfile, {"pH":pH_mat, "redox":redox_mat, "co2":ratio_mat}, oned_as='column')
        
        contour_fig = pylab.figure()
        pylab.hold(True)
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contourf(pH_mat, redox_mat, feasibility_mat, [-1.0, 0.25, 0.75, 2.0])
        CS = plt.contour(pH_mat, redox_mat, ratio_mat, pylab.arange(-4.5, -1.49, 0.5), colors='k')
        plt.clabel(CS, inline=1, fontsize=7, colors='black')
        plt.xlim(min(pH_list), max(pH_list))
        plt.ylim(min(redox_list), max(redox_list))
        plt.xlabel("pH")
        plt.ylabel("$\\log{\\frac{[NADPH]}{[NADP+]}}$")
        plt.title("minimal $\\log{[CO_2]}$ required for feasibility")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=640, height=480)

        #heatmat_fig = plt.figure()
        #plt.pcolor(pH_mat, redox_mat, ratio_mat)
        #plt.colorbar()
        #plt.xlim(min(pH_list), max(pH_list))
        #plt.ylim(min(redox_list), max(redox_list))
        #plt.clim((cid2bounds[11][0], cid2bounds[11][1]+1.0))
        #plt.xlabel("pH")
        #plt.ylabel("$\\log{\\frac{[NADPH]}{[NADP+]}}$")
        #plt.title("minimal $\\log{[CO_2]}$ required for feasibility")
        #self.html_writer.embed_matplotlib_figure(heatmat_fig, width=640, height=480)

    def analyze_standard_conditions(self, key, field_map):
        self.get_conditions(field_map)
        S, rids, fluxes, cids = self.get_reactions(key, field_map)
        self.write_reactions_to_html(S, rids, fluxes, cids, show_cids=False)
        self.thermo.WriteFormationEnergiesToHTML(self.html_writer, cids)


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=['observed_only',
                                   'hatzi_only',
                                   'milo_only',
                                   'milo_merged'],
                          default="milo_merged",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermodynamics_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="The name of the thermodynamics file to load.")
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
    options, _ = MakeOpts().parse_args(sys.argv)
    input_filename = os.path.abspath(options.input_filename)
    output_filename = os.path.abspath(options.output_filename)
    if not os.path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename
    print 'Will write output to %s' % output_filename
    
    thermo_source = options.thermodynamics_source
    assert thermo_source in ('observed_only', 'hatzi_only', 'milo_only', 'milo_merged')
    print 'Will load thermodynamic data from source "%s"' % thermo_source
    
    db_loc = options.db_filename
    print 'Reading from DB %s' % db_loc
    db = SqliteDatabase(db_loc)
    
    public_db_loc = options.kegg_db_filename
    print 'Reading from public DB %s' % public_db_loc
    db_public = SqliteDatabase(public_db_loc)
    
    # dG0 =  -E'*nE*F - R*T*ln(10)*nH*pH
    # Where: 
    #    F  = 0.1 (kJ/mol)/mV
    #    nE - change in e-
    #    nH - change in H+
    #    pH - the conditions in which the E' was measured
    #
    # Ferredoxin  ox/red: E' = -380mV (nE = 1, nH = 0) -> dG0 = 38.0 kJ/mol [1]
    # Ubiqinone   ox/red: E' =  113mV (nE = 2, nH = 2) -> dG0 = -103.2 kJ/mol [1]
    # Menaquinone ox/red: E' =  -74mV (nE = 2, nH = 2) -> dG0 = -65.8 kJ/mol [1]
    #
    # [1] - Thauer 1977
    
    observed_thermo_fname = options.thermodynamics_filename
    print 'Loading observed thermodynamic data from %s' % observed_thermo_fname
    observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
        observed_thermo_fname)
    
    if thermo_source == 'hatzi_only':    
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(db, 'hatzi_thermodynamics')
        thermo.AddPseudoisomer( 139, nH=0,  z=1, nMg=0, dG0=0)      # Ferrodoxin(ox)
        thermo.AddPseudoisomer( 138, nH=0,  z=0, nMg=0, dG0=38.0)   # Ferrodoxin(red)
        thermo.AddPseudoisomer( 399, nH=90, z=0, nMg=0, dG0=0)      # Ubiquinone-10(ox)
        thermo.AddPseudoisomer( 390, nH=92, z=0, nMg=0, dG0=-103.2) # Ubiquinone-10(red)
        thermo.AddPseudoisomer( 828, nH=16, z=0, nMg=0, dG0=0)      # Menaquinone(ox)
        thermo.AddPseudoisomer(5819, nH=18, z=0, nMg=0, dG0=-65.8)  # Menaquinone(red)
        thermo.SetPseudoisomerMap(101, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=0.0)) # THF
        thermo.SetPseudoisomerMap(234, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=-137.5)) # 10-Formyl-THF
        thermo.SetPseudoisomerMap(445, PseudoisomerMap(nH=22, z=0, nMg=0, dG0=65.1)) # 5,10-Methenyl-THF
        thermo.SetPseudoisomerMap(143, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=77.9)) # 5,10-Methylene-THF
        thermo.SetPseudoisomerMap(440, PseudoisomerMap(nH=25, z=0, nMg=0, dG0=32.1)) # 5-Methyl-THF
    elif thermo_source == 'milo_only':
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(
            db, 'gc_pseudoisomers')
    elif thermo_source == 'milo_merged':
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(
            db, 'gc_pseudoisomers')
        thermo.override_data(observed_thermo)
    elif thermo_source == 'observed_only':
        thermo = observed_thermo
    else:
        logging.fatal('Unknown thermodynamic data source.')
    
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

    
