#!/usr/bin/python

import cplex
import logging
import os.path
import pylab
import sys

from pygibbs.stoichiometric_lp import Stoichiometric_LP
from pygibbs.kegg import KeggPathologic
from pygibbs.feasibility import thermodynamic_pathway_analysis
from pygibbs.kegg_utils import write_kegg_pathway
from toolbox.html_writer import HtmlWriter
from toolbox import util
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution

class Pathologic(object):
    
    ALLOWED_THERMODYNAMIC_METHODS = ['none', 'pCr', 'MTDF', 'global', 'localized'] 
    
    def __init__(self, db, public_db, html_writer,
                 thermodynamic_method='global',
                 max_reactions=None,
                 max_solutions=100,
                 maximal_dG=0,
                 update_file='../data/thermodynamics/database_updates.txt'):
        """Initialize the Pathologic object.
        
        Args:
            db: the DB to read group contribution data from.
            html_writer: an HtmlWriter for writing output.
            thermodynamic_method: the analysis methods.
                Options are: "none", "pCr", "MTDF", "global" or "localized"
            max_reactions: the maximum number of reactions to find in a solution.
            max_solutions: the maximum number of solutions to find.
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
        self.thermo = GroupContribution(db=self.db, transformed=False)
        self.thermo.init()
                
        self.flux_relaxtion_factor = None
        self.kegg_patholotic = KeggPathologic()
        self.kegg_patholotic.update_database(update_file, self.html_writer)

    def find_path(self, experiment_name, source=None, target=None):
        """Find a pathway from the source to the target.
        
        Args:    
            experiment_name: a name given to this experiment.
            source: a sparse reaction of the net reactions reactants.
            target: a sparse reaction of the net reactions products.
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
        
        if source is not None:
            exp_html.write('<h2>Source Reaction:</h2>\n')
            exp_html.write_ul(['%d x %s(C%05d)' % (coeff, self.kegg_patholotic.cid2compound[cid].name, cid) for (cid, coeff) in source.iteritems()])
        if target is not None:
            exp_html.write('<h2>Target (biomass) Reaction:</h2>\n')
            exp_html.write_ul(['%d x %s(C%05d)' % (coeff, self.kegg_patholotic.cid2compound[cid].name, cid) for (cid, coeff) in target.iteritems()])

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
        slip.add_stoichiometric_constraints(f, S, compounds, reactions, source, target)
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
        milp.add_stoichiometric_constraints(f, S, compounds, reactions, source, target)
        milp.add_milp_variables()
        if (self.flux_relaxtion_factor != None):
            milp.add_flux_constraint(best_flux * self.flux_relaxtion_factor)
        if (self.max_reactions != None):
            milp.add_reaction_num_constraint(self.max_reactions)
        
        if (self.thermodynamic_method == "pCr"):
            milp.add_dGr_constraints(self.thermo, pCr=True, MTDF=False, maximal_dG=0)
        elif (self.thermodynamic_method == "MTDF"):
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=True, maximal_dG=0)
        elif (self.thermodynamic_method == "global"):
            milp.add_dGr_constraints(self.thermo, pCr=False, MTDF=False, maximal_dG=self.maximal_dG)
        elif (self.thermodynamic_method == "localized"):
            milp.add_localized_dGf_constraints(self.thermo)
        
        for index in range(1, self.max_solutions+1):
            # create the MILP problem to constrain the previous solutions not to reappear again.
            logging.info("Round %03d, solving using MILP" % (milp.solution_index))
            milp.set_objective()
            milp.export("../res/pathologic/%s/%03d_lp.txt" % (experiment_name, milp.solution_index))
            exp_html.write('<b>Solution #%d</b> (<a href="%s/%03d_lp.txt">LP file</a>): '  % (index, experiment_name, index))
            if (not milp.solve()):
                exp_html.write("<b>No solution found</b>")
                logging.info("No more solutions. Quitting!")
                return
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
        Gdot.write('../res/pathologic/%s/%s_graph.svg' % (experiment_name, solution_id), prog='dot', format='svg')
        exp_html.write(' <a href="%s/%s_graph.svg" target="_blank">network</a>' % (experiment_name, solution_id))

        exp_html.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (solution_id))
        exp_html.write('<div id="%s" style="display:none">' % solution_id)

        # write the solution for the concentrations in a table
        if (lp.use_dG_f):
            (cids, concentrations) = lp.get_conc()
            exp_html.write('<p>Compound Concentrations<br>\n')
            exp_html.write('<table border="1">\n')
            exp_html.write('  ' + '<td>%s</td>'*2 % ("KEGG CID", "Concentration [M]") + '\n')
            for c in xrange(len(cids)):
                exp_html.write('<tr><td>C%05d</td><td>%.2g</td></tr>\n' % (cids[c], concentrations[c]))
            exp_html.write('</table></br>\n')
        
        # perform feasibility analysis and write the results
        res = self.margin_analysis(exp_html, sol_reactions, sol_fluxes, experiment_name, solution_id)
        write_kegg_pathway(exp_html, sol_reactions, sol_fluxes)
        exp_html.write('</div>\n')
        
        for optimization in res.keys():
            score = res[optimization][2]
            exp_html.write(", %s = %g" % (optimization, score))
        if not res:
            exp_html.write(", No feasible solutions found")
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

    def margin_analysis(self, exp_html, reactions, fluxes, experiment_name, solution_id):
        cids = [] # I am not using a set() since I want to keep the order of compounds the same as they appear in the reaction
        for r in reactions:
            for cid in r.sparse.keys():
                if (cid not in cids):
                    cids.append(cid)
                
        # convert the list of reactions to a stoichiometric matrix - S
        Nr = len(reactions)
        Nc = len(cids)
        S = pylab.zeros((Nr, Nc))
        rids = []
        for r in range(Nr):
            rids.append(reactions[r].rid)
            for (cid, coeff) in reactions[r].sparse.iteritems():
                c = cids.index(cid)
                S[r, c] = coeff

        return thermodynamic_pathway_analysis(S, rids, fluxes, cids, self.thermo, exp_html)

################################################################################
#                               MAIN                                           #
################################################################################

def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    db = SqliteDatabase('../res/gibbs.sqlite', 'r')
    public_db = SqliteDatabase('../data/public_data.sqlite')
    html_writer = HtmlWriter('../res/pathologic.html')
    update_file = '../data/thermodynamics/database_updates_fermentation.txt'
    pl = Pathologic(db,
                    public_db,
                    html_writer,
                    max_solutions=5,
                    max_reactions=14,
                    thermodynamic_method='MTDF',
                    update_file=update_file)
    pl.thermo.c_range = (1e-6, 1e-2)
    
    # Handy reference
    atp = 2
    adp = 8
    pi = 9
    nad = 3
    nadh = 4
    
    glucose = 31
    g6p = 92
    fbp = 354
    bpg = 236
    g3p = 118
    threepg = 197
    pep = 74
    pyruvate = 22
    succinyl_coa = 91
    acetyl_coa = 24
    lactate = 186
    acetate = 33
    carbon = 19202
    
    source = {glucose:1, adp:3}
    target = {lactate:2, atp:3}
    name = "glu => 2 lactate (3 ATP)"
    pl.thermodynamic_method = 'global'
    
    #source = {glucose:1, adp:2, nad:2}
    #target = {pyruvate:2, atp:2, nadh:2}
    #name = "MTDF Glucose + 2 ADP + 2 NAD => 2 Pyruvate + 2 ATP + 2 NADH (%g - %g, MTDF = %.1f)" % (pl.gc.c_range[0], pl.gc.c_range[1], pl.maximal_dG)
    pl.find_path(name, source, target)

    #source = {}; target = {48:1}
    #name = "=> glyoxylate (%g - %g, MTDF = %.1f)" % (pl.gc.c_range[0], pl.gc.c_range[1], pl.maximal_dG)

    #source = {}; target = {24:1}
    #name = "=> acetyl-CoA (%g - %g, MTDF = %.1f)" % (pl.gc.c_range[0], pl.gc.c_range[1], pl.maximal_dG)

    #source = {}; target = {197:1}
    #for pl.maximal_dG in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 10, 100]:
    #    name = "=> 3PG (%g - %g, MTDF = %.1f)" % (pl.gc.c_range[0], pl.gc.c_range[1], pl.maximal_dG)
    #    logging.info(name)
    #    pl.find_path(name, source, target)
    
    
    #pl.find_path('Glucose to Butanol (Global)', source={31:1}, target={6142:1}, thermo_method="global")

    # TODO: write a clustering algorithm to understand the solutions
    
####################################################################################################

if __name__ == '__main__':
    main()
