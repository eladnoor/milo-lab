#!/usr/bin/python

import cplex
import logging
import os.path

from pygibbs.stoichiometric_lp import Stoichiometric_LP
from pygibbs.kegg import KeggPathologic
from pygibbs.kegg_utils import write_kegg_pathway
from toolbox.html_writer import HtmlWriter
from toolbox import util
from pygibbs.thermodynamic_estimators import LoadAllEstimators

class Pathologic(object):
    
    ALLOWED_THERMODYNAMIC_METHODS = ['none', 'pCr', 'MTDF', 'global', 'localized'] 
    
    def __init__(self, db, public_db, html_writer,
                 thermodynamic_method='global',
                 max_reactions=None,
                 max_solutions=100,
                 maximal_dG=0.0,
                 update_file=None):
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
        assert thermodynamic_method in self.ALLOWED_THERMODYNAMIC_METHODS
        
        cplex.Cplex() # causes CPLEX to print its initialization message
        util._mkdir('../res/pathologic')
        
        self.html_writer = html_writer
        self.thermodynamic_method = thermodynamic_method
        self.max_reactions = max_reactions
        self.max_solutions = max_solutions
        self.maximal_dG = maximal_dG
        
        self.db_public = public_db
        self.db = db
        self.thermo = LoadAllEstimators()['merged']
                
        self.flux_relaxtion_factor = None
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

        exp_html.write('<input type="button" class="button" onclick="return toggleMe(\'__parameters__\')" value="Show Parameters">\n')
        exp_html.write('<div id="__parameters__" style="display:none">')

        exp_html.write('<h2>Conditions:</h2> pH = %g, I = %g, T = %g<br>\n' % (self.thermo.pH, self.thermo.I, self.thermo.T))
        exp_html.write('<h2>Thermodynamic constraints:</h2> ')
        if self.thermodynamic_method == "none":
            exp_html.write("ignore thermodynamics")
        elif self.thermodynamic_method == "pCr":
            exp_html.write("Concentration Range Requirement Analysis, Cmid = %g M" % self.thermo.c_mid)
        elif self.thermodynamic_method == "MTDF":
            exp_html.write("Maximal Chemical Motive Force Analysis, %g M < C < %g M" % self.thermo.c_range)
        elif self.thermodynamic_method == "global":
            exp_html.write("Global constraints, %g M < C < %g M, dG < %.1f" % (self.thermo.c_range[0], self.thermo.c_range[1], self.maximal_dG))
        elif self.thermodynamic_method == "localized":
            exp_html.write("Localized bottlenecks, %g M < C < %g M" % self.thermo.c_range)
        else:
            raise Exception("thermodynamic_method must be one of %s" % self.ALLOWED_THERMODYNAMIC_METHODS)
        exp_html.write('<br>\n')
        
        exp_html.write('<h2>Overall Reaction:</h2>\n')
        exp_html.write(net_reaction.to_hypertext())
        f, S, compounds, reactions = self.kegg_patholotic.get_unique_cids_and_reactions()
        exp_html.write('<h2>%d reactions with %d unique compounds</h2>\n' % (len(reactions), len(compounds)))
        
        exp_html.write('</div><br>\n')
        
        logging.debug("All compounds:")
        for i, compound in enumerate(compounds):
            logging.debug("%05d) C%05d = %s" % (i, compound.cid, compound.name))
        logging.debug("All reactions:")
        for i, reaction in enumerate(reactions):
            logging.debug("%05d) R%05d = %s" % (i, reaction.rid, str(reaction)))

        # Find a solution with a minimal total flux
        logging.info("Preparing the CPLEX object for solving the minimal flux problem")
        exp_html.write('<b>Minimum flux</b>')
        slip = Stoichiometric_LP("Pathologic")
        slip.add_stoichiometric_constraints(f, S, compounds, reactions, net_reaction)
        slip.set_objective()
        slip.export("../res/pathologic/%s/%03d_lp.txt" % (experiment_name, 0))
        exp_html.write(' (<a href="%s/%03d_lp.txt">LP file</a>): ' % (experiment_name, 0))
        logging.info("Solving")
        if not slip.solve():
            exp_html.write("<b>There are no solutions!</b>")
            logging.warning("There are no solutions. Quitting!")
            return
        logging.info("writing solution")
        best_flux = slip.get_total_flux()
        self.write_current_solution(exp_html, slip, experiment_name)

        logging.info("Preparing the CPLEX object for solving the minimal reaction problem using MILP")
        milp = Stoichiometric_LP("Pathologic")
        milp.solution_index = 1
        milp.add_stoichiometric_constraints(f, S, compounds, reactions, net_reaction)
        milp.add_milp_variables()
        if self.flux_relaxtion_factor is not None:
            milp.add_flux_constraint(best_flux * self.flux_relaxtion_factor)
        if self.max_reactions is not None:
            milp.add_reaction_num_constraint(self.max_reactions)
        
        if self.thermodynamic_method == "pCr":
            milp.add_dGr_constraints(self.thermo, pCr=True, MTDF=False, maximal_dG=0)
        elif self.thermodynamic_method == "MTDF":
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=True, maximal_dG=0)
        elif self.thermodynamic_method == "global":
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=False, maximal_dG=self.maximal_dG)
        elif self.thermodynamic_method == "localized":
            milp.add_localized_dGf_constraints(self.thermo)
        
        index = 0
        while (self.max_solutions is None) or (index < self.max_solutions):
            index += 1
            # create the MILP problem to constrain the previous solutions not to reappear again.
            logging.info("Round %03d, solving using MILP" % (milp.solution_index))
            milp.set_objective()
            milp.export("../res/pathologic/%s/%03d_lp.txt" % (experiment_name, milp.solution_index))
            exp_html.write('<b>Solution #%d</b> (<a href="%s/%03d_lp.txt">LP file</a>): '  % (index, experiment_name, index))
            if not milp.solve():
                exp_html.write("<b>No solution found</b>")
                logging.info("No more solutions. Quitting!")
                break
            logging.info("writing solution")
            self.write_current_solution(exp_html, milp, experiment_name)
            milp.ban_current_solution()
        
        exp_html.close()

    def write_current_solution(self, exp_html, lp, experiment_name):
        sol_reactions, sol_fluxes = lp.get_fluxes()
        solution_id = '%03d' % lp.solution_index
        
        exp_html.write('%d reactions, flux = %g, \n' % (len(sol_reactions), lp.get_total_flux()))

        # draw network as a graph and link to it
        Gdot = self.kegg_patholotic.draw_pathway(sol_reactions, sol_fluxes)
        svg_fname = '%s/%s_graph' % (experiment_name, solution_id)
        exp_html.embed_dot_inline(Gdot, width=240, height=320, name=svg_fname)

        if False:
            exp_html.insert_toggle(start_here=True)

            # write the solution for the concentrations in a table
            if lp.use_dG_f:
                cids, concentrations = lp.get_conc()
                exp_html.write('<p>Compound Concentrations<br>\n')
                exp_html.write('<table border="1">\n')
                exp_html.write('  ' + '<td>%s</td>'*2 % ("KEGG CID", "Concentration [M]") + '\n')
                for c in xrange(len(cids)):
                    exp_html.write('<tr><td>C%05d</td><td>%.2g</td></tr>\n' % (cids[c], concentrations[c]))
                exp_html.write('</table></br>\n')
            
            # write the pathway in KEGG format
            write_kegg_pathway(exp_html, sol_reactions, sol_fluxes)
            
            exp_html.div_end()
        exp_html.write('<br>\n')

    def show_Gdot(self, Gdot):
        import gtk
        from toolbox import xdot
    
        win = xdot.DotWindow()
        win.connect('destroy', gtk.main_quit)
        win.set_filter('dot')
        fname = '.dot'
        Gdot.write(fname, format='dot')
        win.open_file(fname)
        gtk.main()

