import sys, pylab, cplex, common
from stoichiometric_lp import Stoichiometric_LP
from kegg import KeggPathologic
from groups import GroupContribution
from feasibility import find_pCr, LinProgNoSolutionException
from html_writer import HtmlWriter

################################################################################
#                               CONSTANTS & DEFAULTS                           #
################################################################################
class Pathologic:
    def __init__(self, update_file='../data/database_updates.txt'):
        cplex.Cplex() # causes CPLEX to print its initialization message

        common._mkdir('../res')
        common._mkdir('../res/pathologic')
        self.LOG_FILE = open('../res/pathologic/pathologic.log', 'w')
        self.UPDATE_FILE = update_file
        self.gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="pathologic", log_file=self.LOG_FILE)
        self.gc.init()
        self.pH = 7
        self.I = 0.1
        self.T = 300
        self.c_range = (1e-6, 1e-2)
        self.c_mid = 1e-4
        self.thermodynamic_method = "margin" # options are: "none", "margin", "global"
        self.max_reactions = None
        self.flux_relaxtion_factor = None
        self.cid2bounds = self.gc.kegg().cid2bounds
        self.cid2dG0_f = self.gc.get_cid2dG0(self.pH, self.I, self.T)
        self.kegg = KeggPathologic(self.LOG_FILE, self.gc.kegg())

    def __del__(self):
        self.LOG_FILE.close()
        
    def find_path(self, experiment_name, source=None, target=None, thermo_method=None, max_reactions=None):
        if (thermo_method != None):
            self.thermodynamic_method = thermo_method
        if (max_reactions != None):
            self.max_reactions = max_reactions
        
        common._mkdir('../res/pathologic/' + experiment_name)
        self.gc.HTML.write('<a href="pathologic/' + experiment_name + '.html">' + experiment_name + '</a><br>\n')
        exp_html = HtmlWriter('../res/pathologic/' + experiment_name + '.html')
        exp_html.write("<p><h1>%s</h1>\n" % experiment_name)

        exp_html.write('<input type="button" class="button" onclick="return toggleMe(\'__parameters__\')" value="Show Parameters">\n')
        exp_html.write('<div id="__parameters__" style="display:none">')

        exp_html.write('<h2>Conditions:</h2> pH = %g, I = %g, T = %g<br>\n' % (self.pH, self.I, self.T))
        exp_html.write('<h2>Thermodynamic constraints:</h2> ')
        if (self.thermodynamic_method == "none"):
            exp_html.write("ignore thermodynamics")
        elif (self.thermodynamic_method == "margin"):
            exp_html.write("Concentration Range Requirement Analysis, Cmid = %g M" % self.c_mid)
        elif (self.thermodynamic_method == "global"):
            exp_html.write("Global constraints, %g M < C < %g M" % self.c_range)
        elif (self.thermodynamic_method == "localized"):
            exp_html.write("Localized bottlenecks, %g M < C < %g M" % self.c_range)
        else:
            raise Exception("thermodynamic_method must be: 'none', 'margin', 'global' or 'localized'")
        exp_html.write('<br>\n')
        
        if (source != None):
            exp_html.write('<h2>Source Reaction:</h2>\n')
            exp_html.write_ul(['%d x %s(C%05d)' % (coeff, self.kegg.cid2compound[cid].name, cid) for (cid, coeff) in source.iteritems()])
        if (target != None):
            exp_html.write('<h2>Target (biomass) Reaction:</h2>\n')
            exp_html.write_ul(['%d x %s(C%05d)' % (coeff, self.kegg.cid2compound[cid].name, cid) for (cid, coeff) in target.iteritems()])
        exp_html.flush()

        self.kegg.update_database(self.UPDATE_FILE, exp_html)
        (f, S, compounds, reactions) = self.kegg.get_unique_cids_and_reactions()
        exp_html.write('<h2>%d reactions with %d unique compounds</h2>\n' % (len(reactions), len(compounds)))
        
        exp_html.write('</div><br>\n')
        
        self.LOG_FILE.write("All compounds:\n")
        for c in range(len(compounds)):
            self.LOG_FILE.write("%05d) C%05d = %s\n" % (c, compounds[c].cid, compounds[c].name))
        self.LOG_FILE.write("All reactions:\n")
        for r in range(len(reactions)):
            self.LOG_FILE.write("%05d) R%05d = %s\n" % (r, reactions[r].rid, str(reactions[r])))

        # Find a solution with a minimal total flux
        sys.stderr.write("Preparing the CPLEX object for solving the minimal flux problem ... ")
        exp_html.write('<b>Minimum flux</b>')
        slip = Stoichiometric_LP("Pathologic", self.LOG_FILE)
        slip.add_stoichiometric_constraints(f, S, compounds, reactions, source, target)
        slip.set_objective()
        slip.export("../res/pathologic/%s/lp_minflux.txt" % experiment_name)
        exp_html.write(' (<a href="%s/lp_minflux.txt">LP file</a>): ' % experiment_name)
        sys.stderr.write("[DONE]\n")
        sys.stderr.write("Solving ... ")
        if (not slip.solve() ):
            exp_html.write("<b>There are no solutions!</b>")
            sys.stderr.write("There are no solutions. Quitting!")
            return
        sys.stderr.write("writing solution ...")
        best_flux = slip.get_total_flux()
        self.write_current_solution(exp_html, slip, experiment_name)
        sys.stderr.write("[DONE]\n")

        sys.stderr.write("Preparing the CPLEX object for solving the minimal reaction problem using MILP ... ")
        milp = Stoichiometric_LP("Pathologic", self.LOG_FILE)
        milp.add_stoichiometric_constraints(f, S, compounds, reactions, source, target)
        milp.add_milp_variables()
        if (self.flux_relaxtion_factor != None):
            milp.add_flux_constraint(best_flux * self.flux_relaxtion_factor)
        if (self.max_reactions != None):
            milp.add_reaction_num_constraint(self.max_reactions)
        
        if (self.thermodynamic_method == "margin"):
            milp.add_dGr_constraints(self.cid2dG0_f)
            milp.add_specific_dGf_constraints(self.cid2bounds)
            milp.add_margin_dGf_constraints(self.c_mid)
        elif (self.thermodynamic_method == "global"):
            milp.add_dGr_constraints(self.cid2dG0_f)
            milp.add_specific_dGf_constraints(self.cid2bounds)
            milp.add_global_dGf_constraints(self.c_range)
        elif (self.thermodynamic_method == "localized"):
            milp.add_localized_dGf_constraints(self.cid2dG0_f, self.cid2bounds, self.c_range)
        
        sys.stderr.write("[DONE]\n")

        while True:
            # create the MILP problem to constrain the previous solutions not to reappear again.
            index = milp.solution_index + 1
            sys.stderr.write("Round %03d, solving using MILP ... " % (index))
            milp.set_objective()
            milp.export("../res/pathologic/%s/lp_%03d.txt" % (experiment_name, index))
            exp_html.write('<b>Solution #%d</b> (<a href="%s/lp_%03d.txt">LP file</a>): '  % (index, experiment_name, index))
            if (not milp.solve()):
                exp_html.write("<b>No solution found</b>")
                sys.stderr.write("No more solutions. Quitting!")
                return
            sys.stderr.write("writing solution ...")
            milp.ban_current_solution()
            self.write_current_solution(exp_html, milp, experiment_name)
            sys.stderr.write("[DONE]\n")
        exp_html.close()

    def write_current_solution(self, exp_html, lp, experiment_name):
        solution = lp.get_fluxes()
        exp_html.write('%d reactions, flux = %g, \n' % (len(solution), lp.get_total_flux()))

        solution_id = 'solution_%03d' % lp.solution_index
        exp_html.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (solution_id))
        exp_html.write('<div id="%s" style="display:none">' % solution_id)

        sol_fluxes = [flux for (r, flux) in solution]
        sol_reactions = [lp.reactions[r] for (r, flux) in solution]
        Gdot = self.kegg.draw_pathway(sol_fluxes, sol_reactions)
        Gdot.write('../res/pathologic/%s/%s_Gdot.svg' % (experiment_name, solution_id), prog='dot', format='svg')
        self.write_kegg_pathway(exp_html, solution, lp.reactions)

        pCr = self.margin_analysis(exp_html, sol_reactions, sol_fluxes, self.cid2dG0_f, experiment_name, solution_id)
            
        exp_html.write('</div>\n')
        exp_html.write(' <a href="%s/%s_Gdot.svg" target="_blank">network</a>' % (experiment_name, solution_id))
        
        if (self.thermodynamic_method == "none"):
            exp_html.write('<br>\n')
        elif (pCr != None):
            exp_html.write(', pCr = %.1f<br>\n' % pCr)
        else:
            exp_html.write(', infeasible<br>\n')

        exp_html.flush()

    def show_Gdot(self, Gdot):
        import xdot, gtk
    
        win = xdot.DotWindow()
        win.connect('destroy', gtk.main_quit)
        win.set_filter('dot')
        fname = '.dot'
        Gdot.write(fname, format='dot')
        win.open_file(fname)
        gtk.main()
    
    def write_kegg_pathway(self, exp_html, solution, reactions):

        def write_reaction(prefix, reaction, flux=1):
            if (flux == 1):
                exp_html.write('%sR%05d&nbsp;&nbsp;%s<br>\n' % (prefix, reaction.rid, str(reaction)))
            else:
                exp_html.write('%sR%05d&nbsp;&nbsp;%s (x%g)<br>\n' % (prefix, reaction.rid, str(reaction), flux))
        
        exp_html.write('<p style="font-family: courier; font-size:10pt">')
        exp_html.write('ENTRY' + '&nbsp;'*7 + 'M-PATHOLOGIC<br>\n')
        exp_html.write('SKIP' + '&nbsp;'*8 + 'FALSE<br>\n')
        exp_html.write('NAME' + '&nbsp;'*8 + 'M-PATHOLOGIC<br>\n')
        exp_html.write('TYPE' + '&nbsp;'*8 + 'MARGIN<br>\n')
        exp_html.write('CONDITIONS' + '&nbsp;'*2 + 'pH=%g,I=%g,T=%g<br>\n' % (self.pH, self.I, self.T))
        exp_html.write('C_MID' + '&nbsp;'*7 + '0.0001<br>\n')
        for i in range(len(solution)):
            (r, flux) = solution[i]
            if (i == 0):
                write_reaction('REACTION' + '&nbsp;'*4, reactions[r], flux)
            else:
                write_reaction('&nbsp;'*12, reactions[r], flux)
        exp_html.write('///<br></p>\n')
        exp_html.flush()

    def margin_analysis(self, exp_html, reactions, fluxes, cid2dG0_f, experiment_name, solution_id):
        
        ###########################
        def dG_to_str(dG):
            if (pylab.isnan(dG)):
                return "N/A"
            else:
                return "%.1f" % dG
        ###########################
        
        cids = [] # I am not using a set() since I want to keep the order of compounds the same as they appear in the reaction
        relevant_reactions = []
        for r in reactions:
            relevant_reactions.append(r)
            for cid in r.sparse.keys():
                if (cid not in cids):
                    cids.append(cid)
                
        # convert the list of reactions to a stoichiometric matrix
        Nr = len(relevant_reactions)
        Nc = len(cids)

        dG0_f = pylab.zeros((Nc, 1))
        ind_nan = []

        # write a table of the compounds and their dG0_f
        exp_html.write('<table border="1">\n')
        exp_html.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f' [kJ/mol]"))
        for c in range(Nc):
            compound = self.kegg.cid2compound[cids[c]]
            cid_str = '<a href="%s">C%05d</a>' % (compound.get_link(), compound.cid)
            
            if cids[c] in cid2dG0_f:
                dG0_f[c, 0] = cid2dG0_f[cids[c]]
            else:
                ind_nan.append(c)
                dG0_f[c, 0] = pylab.nan
            exp_html.write('<tr><td>%s</td><td>%s</td><td>%s</td>\n' % (cid_str, compound.name, dG_to_str(dG0_f[c, 0])))
        exp_html.write('</table><br>\n')

        # write a table of the reactions and their dG0_r
        exp_html.write('<table border="1">\n')
        exp_html.write('  <td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG RID", "Reaction", "flux", "dG0_r' [kJ/mol]"))
        S = pylab.zeros((Nr, Nc))
        dG0_r = pylab.zeros((Nr, 1))
        for r in range(Nr):
            reaction = relevant_reactions[r] 
            rid_str = '<a href="%s">R%05d</a>' % (reaction.get_link(), reaction.rid)
            reaction_str = self.kegg.sparse_to_hypertext(reaction.sparse)
            flux = fluxes[r]

            for (cid, coeff) in reaction.sparse.iteritems():
                c = cids.index(cid)
                S[r, c] = coeff
                dG0_r[r, 0] += coeff*dG0_f[c, 0]
            
            exp_html.write('<tr><td>%s</td><td>%s</td><td>%g</td><td>%s</td>\n' % (rid_str, reaction_str, flux, dG_to_str(dG0_r[r, 0])))
        exp_html.write('</table><br>\n')
          
        bounds = [self.cid2bounds.get(cid, (None, None)) for cid in cids]
        try:
            (dG_f, concentrations, pCr) = find_pCr(S, dG0_f, c_mid=self.c_mid, bounds=bounds)
            dG_r = pylab.dot(S, dG_f)
        except LinProgNoSolutionException:
            exp_html.write('<b>No feasible solution found, cannot calculate the Margin</b>')
            exp_html.write('</p>\n')
            return None

        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 8
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 5
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        pylab.figure()
        pylab.plot(concentrations, range(Nc, 0, -1), '*b')
        pylab.xscale('log')
        pylab.title('Concentrations (pCr = %.1f)' % pCr)
        pylab.ylabel('Compound no.')
        pylab.xlabel('Concentration [M]')

        dG0_f[ind_nan] = dG_f[ind_nan] # since we don't know the dG0_f of some compounds, we will use the solution dG_f instead
        x_min = concentrations.min()/10
        x_max = concentrations.max()*10
        y_min = 0
        y_max = Nc+1
        
        for i in range(Nc):
            pylab.text(concentrations[i, 0]*1.1, Nc-i, self.kegg.cid2compound[cids[i]].name, fontsize=8, rotation=0)
            (b_low, b_up) = bounds[i]
            if (b_low == None):
                b_low = x_min
            if (b_up == None):
                b_up = x_max
            pylab.plot([b_low, b_up], [Nc-i, Nc-i], '-k', linewidth=0.4)
        range_factor = 10**(pCr/2)
        pylab.axvspan(self.c_mid / range_factor, self.c_mid * range_factor, facecolor='g', alpha=0.3)
        pylab.axis([x_min, x_max, y_min, y_max])
        pylab.savefig('../res/pathologic/%s/%s_conc.svg' % (experiment_name, solution_id), format='svg')
        exp_html.embed_svg('%s/%s_conc.svg' % (experiment_name, solution_id), width=800, height=600)

        pylab.figure()
        dG_profile_matrix = pylab.zeros((Nr+1, 3))

        dG_f_mid = dG0_f + common.R * self.T * pylab.log(self.c_mid)
        dG_r_standard = pylab.dot(S, dG0_f)
        dG_r_mid = pylab.dot(S, dG_f_mid)

        dG_profile_matrix[1:,0] = pylab.cumsum(dG_r_standard)
        dG_profile_matrix[1:,1] = pylab.cumsum(dG_r_mid)
        dG_profile_matrix[1:,2] = pylab.cumsum(dG_r)
        
        pylab.plot(range(Nr+1), dG_profile_matrix)
        
        for i in range(Nr):
            pylab.text(i+0.5, pylab.mean(dG_profile_matrix[i:(i+2),0]), "R%05d" % reactions[i].rid,\
                       fontsize=10, rotation=45, horizontalalignment='center', backgroundcolor='white')
        pylab.title('Cumulative Reactions Delta G')
        pylab.ylabel('cumulative dG [kJ/mol]')
        pylab.xlabel('Reaction no.')
        pylab.legend(['Standard [1 M]', 'Physiological [10^(%g) M]' % pylab.log10(self.c_mid), 'Optimized'], fancybox=True, shadow=True, loc="lower left")
        pylab.savefig('../res/pathologic/%s/%s_profile.svg' % (experiment_name, solution_id), format='svg')
        exp_html.embed_svg('%s/%s_profile.svg' % (experiment_name, solution_id), width=800, height=600)
        exp_html.write('</p>\n')
        return pCr

################################################################################
#                               MAIN                                           #
################################################################################

def main():
    pl = Pathologic()
    #pl.find_path('glyoxlyate (CRR, Max_R = 15)', source={}, target={48:1}, thermo_method="margin", max_reactions=15) # => glyoxylate
    
    pl.c_range = (1e-9, 1e-2)
    #pl.find_path('glyoxlyate (Global)', source={}, target={48:1}, thermo_method="global")
    pl.find_path('acetyl-CoA 1nM-10mM', source={}, target={24:1}, thermo_method="global")
    #pl.find_path('3PG 1uM-10mM (no added reactions)', source={}, target={197:1}, thermo_method='global')
    #pl.find_path('3PG (no FDH, no Alanine, no Lactate)', source={}, target={197:1}, thermo_method="global")
    #pl.find_path('Glucose to Butanol (Global)', source={31:1}, target={6142:1}, thermo_method="global")

    
    # TODO: When using "margin" optimization, there is no constraint on the total flux, and that can cause unwanted results (such as futile cycles).
    # one solution for this is to find all futile cycles and remove them in post-processing.
    # Note that these futile cycles are not of the thermodynamically infeasible kind (i.e. they are futile only if you ignore the co-factors).
    
####################################################################################################

if (__name__ == '__main__'):
    main()
