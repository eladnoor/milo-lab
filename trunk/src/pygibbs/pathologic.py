#!/usr/bin/python

import logging
import os.path
import numpy as np

from pygibbs.stoichiometric_lp import Stoichiometric_LP
from pygibbs.kegg import KeggPathologic
from pygibbs.kegg_utils import write_kegg_pathway
from toolbox.html_writer import HtmlWriter
from toolbox import util

class Pathologic(object):
    
    THERMO_METHOD_NONE = 'none'
    THERMO_METHOD_PCR = 'pcr'
    THERMO_METHOD_MTDF = 'mtdf'
    THERMO_METHOD_GLOBAL = 'global'
    THERMO_METHOD_LOCALIZED = 'localized'
    
    ALLOWED_THERMODYNAMIC_METHODS = [THERMO_METHOD_NONE, THERMO_METHOD_PCR,
                                     THERMO_METHOD_MTDF, THERMO_METHOD_GLOBAL,
                                     THERMO_METHOD_LOCALIZED] 
    
    def __init__(self, db, public_db, html_writer,
                 thermo=None,
                 thermodynamic_method='global',
                 max_reactions=None,
                 max_solutions=100,
                 maximal_dG=0.0,
                 update_file=None,
                 output_kegg_file=None):
        """Initialize the Pathologic object.
        
        Args:
            db: the DB to read group contribution data from.
            html_writer: an HtmlWriter for writing output.
            thermodynamic_method: the analysis methods.
                Options are: "none", "pCr", "MTDF", "global" or "localized"
            max_reactions: the maximum number of reactions to find in a solution (use None for unlimited)
            max_solutions: the maximum number of solutions to find (use None for unlimited)
            maximal_dG: the maximum dG allowed.
                Use this to change the thermodynamic constraints to have a different
                MTDF. When set to 0, it is the usual feasibility measure.
            update_file: the file to read for KEGG updates.
        """
        assert thermodynamic_method.lower() in self.ALLOWED_THERMODYNAMIC_METHODS
        
        util._mkdir('../res/pathologic')
        
        self.html_writer = html_writer
        self.thermodynamic_method = thermodynamic_method.lower()
        self.max_reactions = max_reactions
        self.max_solutions = max_solutions
        self.maximal_dG = maximal_dG
        
        self.db_public = public_db
        self.db = db
        self.thermo = thermo
                
        self.kegg_patholotic = KeggPathologic()
        if update_file is not None:
            self.kegg_patholotic.update_database(update_file, self.html_writer)
            
    def add_reaction(self, reaction, weight=1.0):
        self.kegg_patholotic.add_reaction(reaction, weight)
    
    def add_cofactor_reaction(self, reaction):
        self.kegg_patholotic.add_cofactor_reaction(reaction)

    def find_path(self, experiment_name, net_reaction):
        """Find a pathway from the source to the target.
        
        Args:    
            experiment_name: a name given to this experiment.
            net_reaction: a Reaction describing the net reaction for the desired paths
        """
        dirname = os.path.join('../res/pathologic/', experiment_name)
        logging.info('Writing output to: %s' % dirname)
        util._mkdir(dirname)
        
        self.html_writer.write('<a href="pathologic/' + experiment_name + '.html">' + experiment_name + '</a><br>\n')
        exp_html = HtmlWriter('../res/pathologic/' + experiment_name + '.html')
        exp_html.write("<p><h1>%s</h1>\n" % experiment_name)

        exp_html.insert_toggle(div_id="__parameters__", start_here=True,
                               label='Show Parameters')
        exp_html.write('<h2>Thermodynamic constraints:</h2> ')
        if self.thermodynamic_method == Pathologic.THERMO_METHOD_NONE:
            exp_html.write("ignore thermodynamics")
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_PCR:
            exp_html.write("Concentration Range Requirement Analysis, Cmid = %g M" % self.thermo.c_mid)
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_MTDF:
            exp_html.write("Optimized Distributed Bottleneck, %g M < C < %g M" % self.thermo.c_range)
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_GLOBAL:
            exp_html.write("Global constraints, %g M < C < %g M, dG < %.1f" %
                           (self.thermo.c_range[0], self.thermo.c_range[1], self.maximal_dG))
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_LOCALIZED:
            exp_html.write("Localized bottlenecks, %g M < C < %g M" % self.thermo.c_range)
        else:
            raise Exception("thermodynamic_method must be one of %s" % self.ALLOWED_THERMODYNAMIC_METHODS)
        exp_html.write('<br>\n')
        
        f, S, compounds, reactions = self.kegg_patholotic.get_unique_cids_and_reactions()

        exp_html.write('<h2>Conditions:</h2>\n')
        exp_html.write_ul(['pH = %g' % self.thermo.pH,
                           'I = %g' % self.thermo.I,
                           'T = %g' % self.thermo.T,
                           "Max &Delta;<sub>r</sub>G' = %.1f" % self.maximal_dG,
                           'Max no. reactions: %d' % (self.max_reactions or -1),
                           'Max no. solutions: %d' % (self.max_solutions or -1),
                           'Overall Reaction: %s' % net_reaction.to_hypertext(),
                           '%d reactions' % len(reactions),
                           '%d unique compounds' % len(compounds)])

        exp_html.div_end()
        exp_html.write('</br>\n')
        
        logging.debug("All compounds:")
        for i, compound in enumerate(compounds):
            logging.debug("%05d) C%05d = %s" % (i, compound.cid, compound.name))
        logging.debug("All reactions:")
        for i, reaction in enumerate(reactions):
            logging.debug("%05d) R%05d = %s" % (i, reaction.rid, str(reaction)))

        output_kegg_file = open(dirname + '/kegg_pathway.txt', 'w')
        exp_html.write('<a href="%s/kegg_pathway.txt">All solutions in KEGG format</a></br>\n'
                       % experiment_name)
        
        # Find a solution with a minimal total flux
        logging.info("Preparing LP solver for the minimal total flux problem")
        exp_html.write('<b>Minimum flux</b>')
        slip = Stoichiometric_LP("Pathologic")
        slip.add_stoichiometric_constraints(f, S, compounds, reactions, net_reaction)
        slip.export("../res/pathologic/%s/%03d_lp.txt" % (experiment_name, 0))
        exp_html.write(' (<a href="%s/%03d_lp.txt">LP file</a>): ' % (experiment_name, 0))
        logging.info("Solving")
        if not slip.solve():
            exp_html.write("<b>There are no solutions!</b>")
            logging.warning("There are no solutions. Quitting!")
            return
        logging.info("writing solution")
        self.write_current_solution(exp_html, slip, experiment_name)

        logging.info("Preparing MILP solver for the minimal no. reaction problem")
        milp = Stoichiometric_LP("Pathologic")
        milp.solution_index = 1
        milp.add_stoichiometric_constraints(f, S, compounds, reactions, net_reaction)
        milp.add_milp_variables()
        if self.max_reactions is not None:
            milp.add_reaction_num_constraint(self.max_reactions)
       
        if self.thermodynamic_method == Pathologic.THERMO_METHOD_PCR:
            milp.add_dGr_constraints(self.thermo, pCr=True, MTDF=False, maximal_dG=0)
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_MTDF:
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=True, maximal_dG=0)
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_GLOBAL:
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=False, maximal_dG=self.maximal_dG)
        elif self.thermodynamic_method == Pathologic.THERMO_METHOD_LOCALIZED:
            milp.add_localized_dGf_constraints(self.thermo)
        
        index = 0
        while (self.max_solutions is None) or (index < self.max_solutions):
            index += 1
            # create the MILP problem to constrain the previous solutions not to reappear again.
            logging.info("Round %03d, solving using MILP" % (milp.solution_index))
            milp.export("../res/pathologic/%s/%03d_lp.txt" % (experiment_name, milp.solution_index))
            exp_html.write('<b>Solution #%d</b> (<a href="%s/%03d_lp.txt">LP file</a>): '  % (index, experiment_name, index))
            if not milp.solve():
                exp_html.write("<b>No solution found</b>")
                logging.info("No more solutions. Quitting!")
                break
            logging.info("writing solution")
            self.write_current_solution(exp_html, milp, experiment_name,
                                        output_kegg_file)
            milp.ban_current_solution()
        
        output_kegg_file.close()
        exp_html.close()

    def write_current_solution(self, exp_html, lp, experiment_name,
                               output_kegg_file=None):
        solution = lp.get_active_reaction_data()
        solution_id = '%03d' % lp.solution_index
        
        exp_html.write('%d reactions, flux = %g, \n' %
                       (len(solution.reactions),
                        float(solution.fluxes.sum(1))))

        # draw network as a graph and link to it
        Gdot = self.kegg_patholotic.draw_pathway(solution.reactions,
                                                 list(solution.fluxes.flat))
        svg_fname = '%s/%s_graph' % (experiment_name, solution_id)
        exp_html.embed_dot_inline(Gdot, width=240, height=320, name=svg_fname)

        # write the solution for the concentrations in a table
        if solution.concentrations is not None:
            exp_html.insert_toggle(start_here=True)
            
            exp_html.write('Compound Concentrations<br>\n')
            rowdicts = []
            for c, compound in enumerate(solution.compounds):
                rowdict = {'KEGG ID':'<a href="%s">C%05d</a>' %
                           (compound.get_link(), compound.cid),
                           'Compound':compound.name}
                if np.isfinite(solution.concentrations[0, c]):
                    rowdict["dG'0 [kJ/mol]"] = '%.1f' % solution.dG0_f[0, c]
                    rowdict['Conc. [M]'] = '%.2g' % solution.concentrations[0, c]
                else:
                    rowdict["dG'0 [kJ/mol]"] = 'N/A'
                    rowdict['Conc. [M]'] = 'N/A'
                rowdict["dG' [kJ/mol]"] = '%.1f' % solution.dG_f[0, c]
                rowdicts.append(rowdict)
            exp_html.write_table(rowdicts,
                headers=['KEGG ID', 'Compound', "dG'0 [kJ/mol]", "dG' [kJ/mol]",
                         'Conc. [M]'])
            
            exp_html.write('Reaction Gibbs energies<br>\n')
            rowdicts = []
            for r, reaction in enumerate(solution.reactions):
                rowdict = {'KEGG ID':'<a href="%s">R%05d</a>' %
                           (reaction.get_link(), reaction.rid),
                           'Reaction':reaction.to_hypertext(show_cids=False)}
                if np.isfinite(solution.dG_r[0, r]):
                    rowdict["dG' [kJ/mol]"] = '%.1f' % solution.dG_r[0, r]
                else:
                    rowdict["dG' [kJ/mol]"] = 'N/A'

                if np.isfinite(solution.dG0_r[0, r]):
                    rowdict["dG'0 [kJ/mol]"] = '%.1f' % solution.dG0_r[0, r]
                else:
                    rowdict["dG'0 [kJ/mol]"] = 'N/A'

                rowdicts.append(rowdict)
            exp_html.write_table(rowdicts, 
                headers=['KEGG ID', 'Reaction', "dG'0 [kJ/mol]", "dG' [kJ/mol]"])

            exp_html.div_end()

        # write the pathway in KEGG format
        if output_kegg_file is not None:
            write_kegg_pathway(output_kegg_file,
                               entry=experiment_name + ' ' + solution_id,
                               reactions=solution.reactions,
                               fluxes=list(solution.fluxes.flat))
            
        exp_html.write('<br>\n')

