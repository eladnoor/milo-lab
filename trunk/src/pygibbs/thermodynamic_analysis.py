#!/usr/bin/python

import logging
import matplotlib
import os
import re
import sys 
import numpy as np

from copy import deepcopy
from optparse import OptionParser
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.kegg import Kegg
from pygibbs.pathway import PathwayData
from pygibbs.thermodynamic_constants import transform, RedoxCarriers
from pygibbs.thermodynamic_constants import default_T, R, F
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from toolbox.util import _mkdir
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pygibbs.pathway_modelling import KeggPathway, DeltaGNormalization
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.compound_abundance import CompoundAbundance


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

            function_dict = {'MTDF':self.analyze_mtdf,
                             'MTDF2D':self.analyze_mtdf_2d,
                             'REDOX':self.analyze_redox,
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
        self.thermo.SetConditionsToDefault()
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        
        self.html_writer.write("<li>Reactions:</br><ul>\n")
        
        for r in range(S.shape[0]):
            self.html_writer.write('<li><a href=' + self.kegg.rid2link(rids[r]) + 
                                   '>%s ' % self.kegg.rid2string(rids[r]) + '</a>')
            self.html_writer.write("[x%g, &Delta;<sub>r</sub>G'&deg; = %.1f] : " % (fluxes[r], dG0_r[r, 0]))
            self.html_writer.write(self.kegg.vector_to_hypertext(S[r, :].flat, cids, show_cids=show_cids))
            self.html_writer.write('</li>\n')
        
        v_total = np.dot(np.matrix(fluxes), S).flat
        dG0_total = np.dot(np.matrix(fluxes), dG0_r)[0,0]
        self.html_writer.write('<li><b>Total </b>')
        self.html_writer.write('[&Delta;<sub>r</sub>G&deg; = %.1f kJ/mol] : \n' % dG0_total)
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
                d["&Delta;<sub>f</sub>G'&deg;"] = 'N/A'
            else:
                d["&Delta;<sub>f</sub>G'&deg;"] = '%.1f' % dG0_f_prime[i, 0]
            dict_list.append(d)
        self.html_writer.write_table(dict_list, 
                headers=['CID', 'name', "&Delta;<sub>f</sub>G'&deg;"])

    def get_conditions(self, pathway_data):
        self.thermo.SetConditions(pH=pathway_data.pH,
                                  I=pathway_data.I,
                                  T=pathway_data.T,
                                  pMg=pathway_data.pMg)
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

    def analyze_mtdf(self, key, pathway_data):
        self.get_conditions(pathway_data)
        cid2bounds = self.get_bounds(key, pathway_data)
        #self.write_bounds_to_html(cid2bounds, self.thermo.c_range)
        S, rids, fluxes, cids = self.get_reactions(key, pathway_data)
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        fluxes = np.matrix(fluxes)
        
        keggpath = KeggPathway(S, rids, fluxes, cids, None, dG0_r,
                               cid2bounds=cid2bounds, c_range=self.thermo.c_range)
        keggpath.normalization = DeltaGNormalization.SIGN_FLUX
        #try:
        _,  mtdf = keggpath.FindMtdf()
        ln_conc, total_dG_prime = keggpath.GetTotalReactionEnergy(mtdf, maximize=True)
            
        #except UnsolvableConvexProblemException as e:
        #    self.html_writer.write("<b>WARNING: cannot calculate MTDF "
        #                           "because %s:</b></br>\n" %
        #                           str(e))
        #    problem_str = str(e.problem).replace('\n', '</br>\n')
        #    self.html_writer.write("%s" % problem_str)
        #    return
        
        odfe = 100 * np.tanh(mtdf / (2*R*self.thermo.T))
        concentrations = np.exp(ln_conc)
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, name=key+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        concentration_fig = keggpath.PlotConcentrations(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=concentration_fig)
        self.html_writer.embed_matplotlib_figure(concentration_fig, name=key+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)

        self.write_metabolic_graph(key+"_grph", S, rids, cids)
        
        #self.write_formation_energies_to_html(cids)
        #dG_r_prime = keggpath.CalculateReactionEnergiesUsingConcentrations(concentrations)
        #return "ODFE = %.1f%%, Total &Delta;<sub>r</sub>G' = %.1f [min = %.1f, max = %.1f] kJ/mol" % \
        #    (odfe, float(np.dot(dG_r_prime.T, fluxes)), 
        #     min_total_dG_prime, max_total_dG_prime)
        
        average_dG_prime = total_dG_prime/np.sum(fluxes)
        average_dfe = 100 * np.tanh(-average_dG_prime / (2*R*self.thermo.T))
        
        print ','.join([key, '%.1f' % mtdf, '%.1f' % -average_dG_prime, 
                        '%.1f' % odfe, '%.1f' % average_dfe, 
                        '%.1f' % total_dG_prime, '%g' % np.sum(fluxes)])
        return "MTDF = %.1f (avg. = %.1f) kJ/mol, ODFE = %.1f%% (avg. = %.1f%%), Total &Delta;<sub>r</sub>G' = %.1f kJ/mol, no. steps = %g" %\
            (mtdf, -average_dG_prime, odfe, average_dfe, total_dG_prime, np.sum(fluxes))
        
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
        self.html_writer.embed_matplotlib_figure(protonation_fig)
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for (nH, z, dG0) in pmatrix:
            self.html_writer.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        self.html_writer.write('</table>')

    def analyze_redox(self, key, pathway_data):
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
        E_list = pathway_data.redox_values or np.arange(-0.500, -0.249999, 0.01)
        
        pH_mat = np.zeros((len(pH_list), len(E_list)))
        E_mat = np.zeros((len(pH_list), len(E_list)))
        ratio_mat = np.zeros((len(pH_list), len(E_list)))
        
        cid_to_minimize = 11 # CO2
        max_allowed_concentration = 1.0 # 1M
        cid2bounds[cid_to_minimize] = (1e-50, max_allowed_concentration)
        
        for i, pH in enumerate(pH_list):
            self.thermo.SetConditions(pH=pH)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, redox in enumerate(E_list):
                pH_mat[i, j] = pH
                E_mat[i, j] = redox
                
                E0 = -0.32 # for NADP+/NADPH (in V)
                r = np.exp(-2*F/(R*default_T) * (redox - E0))
                cid2bounds[6] = (1e-5, 1e-5) # NADP+
                cid2bounds[5] = (1e-5*r, 1e-5*r) # NADPH
                logging.debug("E = %g, ratio = %.1g" % (redox, r))
                
                keggpath = KeggPathway(S, rids, fluxes, cids,
                                       formation_energies=dG0_f,
                                       cid2bounds=cid2bounds,
                                       c_range=self.thermo.c_range)
                #try:
                # minimize CO2 concentration
                _ln_conc, min_ln_co2_conc = keggpath.FindMinimalFeasibleConcentration(cid_to_minimize)
                ratio_mat[i, j] = min_ln_co2_conc / np.log(10) # change to base 10
                #except UnsolvableConvexProblemException as e:
                #    logging.debug(str(e))
                #    ratio_mat[i, j] = np.log10(max_allowed_concentration) + \
                #                      1.0 # i.e. 10 times higher than the upper bound
        
        field_map = pathway_data.field_map  
        matfile = field_map.GetStringField("MATFILE", "")
        if matfile:
            scipy.io.savemat(matfile, {"pH":pH_mat, "redox":E_mat, 
                "co2":ratio_mat, "S":S, "rids":np.array(rids),
                "fluxes":np.array(fluxes), "cids":np.array(cids),
                "dG0_f":dG0_f}, oned_as='column')
        
        fig = plt.figure()
        fig.hold(True)
        #CS = plt.contour(pH_mat, E_mat, ratio_mat, np.arange(-4.5, 0.01, 0.5), colors='k')
        #plt.clabel(CS, inline=1, fontsize=7, colors='black')
        collection = plt.pcolor(pH_mat, E_mat, ratio_mat, figure=fig)
        plt.colorbar(mappable=collection)
        plt.xlim(min(pH_list), max(pH_list))
        plt.ylim(min(E_list), max(E_list))
        plt.xlabel("pH")
        plt.ylabel("$E^'$ (V)")
        plt.title("minimal $\\log([CO_2]/1M)$ required for feasibility")
        self.html_writer.embed_matplotlib_figure(fig)

    def analyze_mtdf_2d(self, key, pathway_data, contour=False):
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

        E = -0.32
        dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
        dG0_f[cids.index(28)] = 0
        dG0_f[cids.index(30)] = -E * F
        keggpath = KeggPathway(S, rids, fluxes, cids,
                               formation_energies=dG0_f,
                               cid2bounds=cid2bounds,
                               c_range=self.thermo.c_range)
        #try:
        ln_conc, mtdf = keggpath.FindMtdf()
        #except UnsolvableConvexProblemException as e:
        #    self.html_writer.write("<b>WARNING: cannot calculate MTDF "
        #                           "because %s:</b></br>\n" %
        #                           str(e))

        odfe = 100 * np.tanh(mtdf / (2*R*self.thermo.T))
        concentrations = np.exp(ln_conc)
        
        profile_fig = keggpath.PlotProfile(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=profile_fig)
        self.html_writer.embed_matplotlib_figure(profile_fig, name=key+"_prfl")
        keggpath.WriteProfileToHtmlTable(self.html_writer, concentrations)
        
        concentration_fig = keggpath.PlotConcentrations(concentrations)
        plt.title('ODFE = %.1f%%' % odfe, figure=concentration_fig)
        self.html_writer.embed_matplotlib_figure(concentration_fig, name=key+"_conc")
        keggpath.WriteConcentrationsToHtmlTable(self.html_writer, concentrations)
        
        self.html_writer.write('</br>\n')
        
        min_pH, max_pH = 5.0, 9.0
        min_E, max_E = -0.5, -0.1
        steps = np.r_[0:1:21j] # array of 21 numbers evenly spaced from 0 to 1

        pH_list = pathway_data.pH_values or (min_pH + steps*(max_pH-min_pH)) 
        E_list = pathway_data.redox_values or (min_E + steps*(max_E-min_E)) 
        
        pH_mat = np.zeros((len(pH_list), len(E_list)))
        E_mat = np.zeros((len(pH_list), len(E_list)))
        odfe_mat = np.zeros((len(pH_list), len(E_list)))
        
        for i, pH in enumerate(pH_list):
            self.thermo.SetConditions(pH=pH)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            for j, E in enumerate(E_list):
                pH_mat[i, j] = pH
                E_mat[i, j] = E
                
                # override the formation energy of the reduced species to match
                # the value of E.
                dG0_f[cids.index(28)] = 0
                dG0_f[cids.index(30)] = -E * F
                keggpath = KeggPathway(S, rids, fluxes, cids,
                                       formation_energies=dG0_f,
                                       cid2bounds=cid2bounds,
                                       c_range=self.thermo.c_range)
                #try:
                _, mtdf = keggpath.FindMtdf()
                odfe_mat[i, j] = 100 * np.tanh(mtdf / (2*R*self.thermo.T))
                    
                #except UnsolvableConvexProblemException as e:
                #    logging.debug(str(e))
                #    odfe_mat[i, j] = np.NaN
        
        if contour:
            fig = plt.figure()
            CS = plt.contour(pH_mat, E_mat, odfe_mat, np.r_[0:100:6j], 
                             colors='k', figure=fig)
            plt.clabel(CS, inline=1, fontsize=7, colors='black', figure=fig)
        else:
            fig = plt.figure(figsize=(8.0, 6.0)) # make more room for colorbar
            cdict = {'red': ((0.0, 1.0, 1.0),
                            (1.0, 0.0, 0.0)),
                    'green': ((0.0, 0.0, 0.0),
                              (1.0, 1.0, 1.0)),
                    'blue': ((0.0, 0.0, 0.0),
                             (1.0, 0.0, 0.0))}
            cmap = colors.LinearSegmentedColormap('red_green', cdict, 256)
            plt.contourf(pH_mat, E_mat, odfe_mat,
                         levels=[0, 20, 40, 60, 80, 100], 
                         cmap=cmap, figure=fig)
            plt.colorbar()
        fig.axes[0].set_xlim(min(pH_list), max(pH_list))
        fig.axes[0].set_ylim(min(E_list), max(E_list))
        plt.xlabel("pH", figure=fig)
        plt.ylabel("$E^'$ (V)", figure=fig)
        plt.title("ODFE (in %)", figure=fig)
        
        redox = RedoxCarriers()
        for name in ['NADP', 'ferredoxin', 'FAD']:
            rc = redox[name]
            E_min = rc.E_prime - 0.5*R*default_T*np.log(10)/(F*rc.delta_e)
            E_max = rc.E_prime + 0.5*R*default_T*np.log(10)/(F*rc.delta_e)
            if min(E_list) > E_max or max(E_list) < E_min:
                continue
            plt.axhspan(ymin=E_min, ymax=E_max, facecolor='w', 
                        alpha=0.5, figure=fig)
            #plt.text(max(pH_list)-10, rc.E_prime, name,
            #         alpha=1.0, ha='right', va='center')
        
        self.html_writer.embed_matplotlib_figure(fig, name=key+"_cntr")

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

    
