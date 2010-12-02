import pylab
import cvxopt.solvers
import csv, sys
from pygibbs.thermodynamics import MissingCompoundFormationEnergy
from matplotlib.font_manager import FontProperties
from toolbox.html_writer import HtmlWriter
from pygibbs import kegg

try:
    import cplex
    IsCplexInstalled = True
except ImportError:
    IsCplexInstalled = False
    
R = 8.31e-3 # gas constant (kJ/K mol)
T = 300 # temperature (K)

class LinProgNoSolutionException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def linprog(f, A, b, lb=[], ub=[], log_stream=None):
    """
        The constraints are:
            A*x <= b
            lb <= x <= ub
        
        The optimization rule is:
            minimize f*x
        
        In matrix A:
            rows are reactions (indexed by 'r') - total of Nr
            columns are compounds (indexed by 'c') - total of Nc
        
        All other parameters (f, b, lb, ub) should be column vectors, with the following sizes:
            f  - Nc x 1
            b  - Nr x 1
            lb - list of pairs (column_index, lower_bound)
            ub - list of pairs (column_index, upper_bound)
    """
    (Nr, Nc) = A.shape
    if (f.shape[0] != Nc or f.shape[1] != 1):
        raise Exception("linprog: 'f' must be a column vector whose length matches the number of columns in 'A'")
    if (b.shape[0] != Nr or b.shape[1] != 1):
        raise Exception("linprog: 'b' must be a column vector whose length matches the number of rows in 'A'")
    
    if False: # currently CPLEX does not work (claims every problem is infeasible)
        for (c, bound) in lb:
            row_lower = pylab.zeros(Nc)
            row_lower[c] = -1
            A = pylab.vstack([A, row_lower])
            b = pylab.vstack([b, [-bound]])
        for (c, bound) in ub:
            row_upper = pylab.zeros(Nc)
            row_upper[c] = 1
            A = pylab.vstack([A, row_upper])
            b = pylab.vstack([b, [bound]])
        
        # make sure the values are floating point (not integer)
        # and convert the matrices to cvxopt format
        c = cvxopt.matrix(f * 1.0)
        G = cvxopt.matrix(A * 1.0)
        h = cvxopt.matrix(b * 1.0)
        cvxopt.solvers.options['show_progress'] = False
        cvxopt.solvers.options['LPX_K_MSGLEV'] = 0
        sol = cvxopt.solvers.lp(c, G, h, solver='glpk')
        if (not sol['x']):
            return None
        else:
            return pylab.matrix(sol['x'])
    else:
        cpl = cplex.Cplex()
        if (log_stream != None):
            cpl.set_log_stream(log_stream)
        else:
            cpl.set_log_stream(None)
        cpl.set_results_stream(None)
        #cpl.set_warning_stream(None)
        
        cpl.set_problem_name('LP')
        cpl.variables.add(names=["c%d" % c for c in range(Nc)], lb=[-1e6]*Nc, ub=[1e6]*Nc)
        cpl.variables.set_lower_bounds(lb)
        cpl.variables.set_upper_bounds(ub)
    
        cpl.objective.set_linear([(c, f[c, 0]) for c in range(Nc)])
        
        for r in range(Nr):
            cpl.linear_constraints.add(senses='L', names=["r%d" % r])
            for c in range(Nc):
                cpl.linear_constraints.set_coefficients(r, c, A[r, c])
        cpl.linear_constraints.set_rhs([(r, b[r, 0]) for r in range(Nr)])
        cpl.solve()
        if (cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.optimal):
            return None
        else:
            return pylab.matrix(cpl.solution.get_values()).T

def create_cplex(S, dG0_f, log_stream=None):
    cpl = cplex.Cplex()
    if (log_stream != None):
        cpl.set_log_stream(log_stream)
    else:
        cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    #cpl.set_warning_stream(None)
    
    Nr, Nc = S.shape
    cpl.set_problem_name('LP')
    cpl.variables.add(names=["c%d" % c for c in range(Nc)], lb=[-1e6]*Nc, ub=[1e6]*Nc)

    for r in xrange(Nr):
        c_list = pylab.find(S[r, :])
        cpl.linear_constraints.add(senses='L', names=["r%d" % r], rhs=[0])
        for c in c_list:
            cpl.linear_constraints.set_coefficients("r%d" % r, "c%d" % c, S[r, c])

    return cpl

def add_thermodynamic_constraints(cpl, dG0_f, c_range=(1e-6, 1e-2), bounds=None):   
    """
        For any compound that does not have an explicit bound set by the 'bounds' argument,
        create a bound using the 'margin' variables (the last to columns of A).
    """
    
    Nc = dG0_f.shape[0]

    if (bounds != None and len(bounds) != Nc):
        raise Exception("The concentration bounds list must be the same length as the number of compounds")
    
    for c in xrange(Nc):
        if (pylab.isnan(dG0_f[c, 0])):
            continue # unknown dG0_f - cannot bound this compound's concentration at all

        if (bounds == None or bounds[c][0] == None):
            # lower bound: dG0_f + R*T*ln(Cmin) <= x_i
            cpl.variables.set_lower_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(c_range[0]))
        else:
            # this compound has a specific lower bound on its activity
            cpl.variables.set_lower_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(bounds[c][0]))

        if (bounds == None or bounds[c][1] == None):
            # upper bound: x_i <= dG0_f + R*T*ln(Cmax)
            cpl.variables.set_upper_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(c_range[1]))
        else:
            # this compound has a specific upper bound on its activity
            cpl.variables.set_upper_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(bounds[c][1]))

def find_mcmf(S, dG0_f, c_range=(1e-6, 1e-2), bounds=None, log_stream=None):
    """
        Find a distribution of concentration that will satisfy the 'relaxed' thermodynamic constraints.
        The 'relaxation' means that there is a slack variable 'B' where all dG_r are constrained to be < B.
        Note that B can also be negative, which will happen when the pathway is feasible.
        MCMF (Maximal Chemical Motive Force) is defined as the minimal B, note that it is a function of the concentration bounds.
    """
    (Nr, Nc) = S.shape
    
    # compute right hand-side vector - r,
    # i.e. the deltaG0' of the reactions divided by -RT
    if (S.shape[1] != dG0_f.shape[0]):
        raise Exception("The S matrix has %d columns, while the dG0_f vector has %d" % (S.shape[1], dG0_f.shape[0]))
    
    cpl = create_cplex(S, dG0_f, log_stream)
    add_thermodynamic_constraints(cpl, dG0_f, c_range, bounds)

    # Define the MCMF variable and use it relax the thermodynamic constraints on each reaction
    cpl.variables.add(names=["B"], lb=[-1e6], ub=[1e6])
    for r in xrange(Nr):
        cpl.linear_constraints.set_coefficients("r%d" % r, "B", -1)

    # Set 'B' as the objective
    cpl.objective.set_linear([("B", 1)])
    
    #cpl.write("../res/test_MCMF.lp", "lp")
    cpl.solve()
    if (cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.optimal):
        raise LinProgNoSolutionException("")

    dG_f = pylab.matrix(cpl.solution.get_values(["c%d" % c for c in xrange(Nc)])).T
    concentrations = pylab.exp((dG_f-dG0_f)/(R*T))
    MCMF = cpl.solution.get_values(["B"])[0]

    return (dG_f, concentrations, MCMF)

def pC_to_range(pC, c_mid=1e-3, ratio=3.0):
    return (c_mid * 10 ** (-ratio/(ratio+1.0) * pC), c_mid * 10 ** (1.0/(ratio+1.0) * pC))

def generate_constraints_pCr(cpl, dG0_f, c_mid, ratio, bounds=None):   
    """
        For any compound that does not have an explicit bound set by the 'bounds' argument,
        create a bound using the 'margin' variables (the last to columns of A).
        
        The constraints will eventually look like this:
        
        ln(c_mid) - r/(1+r) * pC <= ln(c) <= ln(c_mid) + 1/(1+r) * pC 
    """
    
    Nc = dG0_f.shape[0]

    if (bounds != None and len(bounds) != Nc):
        raise Exception("The concentration bounds list must be the same length as the number of compounds")
    
    cpl.variables.add(names=['pC'], lb=[0], ub=[1e6])
    for c in xrange(Nc):
        if (pylab.isnan(dG0_f[c, 0])):
            continue # unknown dG0_f - cannot bound this compound's concentration at all

        if (bounds == None or bounds[c][0] == None):
            # lower bound: x_i + r/(1+r) * R*T*ln(10)*pC >= dG0_f + R*T*ln(Cmid) 
            cpl.linear_constraints.add(senses='G', names=['c%d_lower' % c], rhs=[dG0_f[c, 0] + R*T*pylab.log(c_mid)])
            cpl.linear_constraints.set_coefficients('c%d_lower' % c, 'c%d' % c, 1)
            cpl.linear_constraints.set_coefficients('c%d_lower' % c, 'pC', R*T*pylab.log(10) * ratio / (ratio + 1.0))
        else:
            # this compound has a specific lower bound on its activity
            cpl.variables.set_lower_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(bounds[c][0]))

        row_upper = pylab.zeros(Nc + 1)
        row_upper[c] = 1
        if (bounds == None or bounds[c][1] == None):
            # upper bound: x_i - 1/(1+r) * R*T*ln(10)*pC <= dG0_f + R*T*ln(Cmid)
            cpl.linear_constraints.add(senses='L', names=['c%d_upper' % c], rhs=[dG0_f[c, 0] + R*T*pylab.log(c_mid)])
            cpl.linear_constraints.set_coefficients('c%d_upper' % c, 'c%d' % c, 1)
            cpl.linear_constraints.set_coefficients('c%d_upper' % c, 'pC', -R*T*pylab.log(10) / (ratio + 1.0))
        else:
            # this compound has a specific upper bound on its activity
            cpl.variables.set_upper_bounds('c%d' % c, dG0_f[c, 0] + R*T*pylab.log(bounds[c][1]))

    # Set 'pC' as the objective
    cpl.objective.set_linear([("pC", 1)])

def find_pCr(S, dG0_f, c_mid=1e-3, ratio=3.0, bounds=None, log_stream=None):
    """
        Compute the feasibility of a given set of reactions
    
        input: s1 = stoichiometric matrix (reactions x compounds)
               dG = deltaG0'-formation values for all compounds (in kJ/mol) (1 x compounds)
               c_mid = the default concentration
               ratio = the ratio between the distance of the upper bound from c_mid and the lower bound from c_mid (in logarithmic scale)
        
        output: (concentrations, margin)
    """
    Nc = S.shape[1]
    if (Nc != dG0_f.shape[0]):
        raise Exception("The S matrix has %d columns, while the dG0_f vector has %d" % (Nc, dG0_f.shape[0]))

    cpl = create_cplex(S, dG0_f, log_stream)
    generate_constraints_pCr(cpl, dG0_f, c_mid, ratio, bounds)

    #cpl.write("../res/test_PCR.lp", "lp")
    cpl.solve()
    if (cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.optimal):
        raise LinProgNoSolutionException("")
    dG_f = pylab.matrix(cpl.solution.get_values(["c%d" % c for c in xrange(Nc)])).T
    concentrations = pylab.exp((dG_f-dG0_f)/(R*T))
    pCr = cpl.solution.get_values(["pC"])[0]

    return (dG_f, concentrations, pCr)

def find_unfeasible_concentrations(S, dG0_f, c_range, c_mid=1e-4, bounds=None, log_stream=None):
    """ 
        Almost the same as find_pCr, but adds a global restriction on the concentrations (for compounds
        that don't have specific bounds in 'bounds').
        After the solution which optimizes the pCr is found, any concentration which does not confer
        to the limits of c_range will be truncated to the closes allowed concentration.
        If at least one concentration needs to be adjusted, then pCr looses its meaning
        and therefore is returned with the value None.
    """
    (dG_f, concentrations, pCr) = find_pCr(S, dG0_f, c_mid, bounds, log_stream)

    for c in xrange(dG0_f.shape[0]):
        if (pylab.isnan(dG0_f[c, 0])):
            continue # unknown dG0_f - therefore the concentration of this compounds is meaningless

        if ((bounds == None or bounds[c][0] == None) and concentrations[c, 0] < c_range[0]):
            concentrations[c, 0] = c_range[0]
            dG_f[c, 0] = dG0_f[c, 0] + R * T * c_range[0]
            pCr = None
        elif ((bounds == None or bounds[c][1] == None) and concentrations[c, 0] > c_range[1]):
            concentrations[c, 0] = c_range[1]
            dG_f[c, 0] = dG0_f[c, 0] + R * T * c_range[1]
            pCr = None

    return (dG_f, concentrations, pCr)

def thermodynamic_pathway_analysis(S, rids, fluxes, cids, thermodynamics, kegg, html_writer):
    (Nr, Nc) = S.shape

    # adjust the directions of the reactions in S to fit the fluxes
    fluxes = [abs(f) for f in fluxes]
    
    kegg.write_reactions_to_html(html_writer, S, rids, fluxes, cids, show_cids=False)

    # calculate the dG0_f of each compound
    dG0_f = pylab.zeros((Nc, 1))
    ind_nan = []
    html_writer.write('<table border="1">\n')
    html_writer.write('  ' + '<td>%s</td>'*5 % ("KEGG CID", "Compound Name", "dG0_f [kJ/mol]", "nH", "z") + '\n')
    for c in range(Nc):
        cid = cids[c]
        name = kegg.cid2name(cid)
        try:
            for (nH, z, dG0) in  thermodynamics.cid2pmatrix(cid):
                html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % \
                                  (kegg.cid2link(cid), cid, name, dG0, nH, z))
            dG0_f[c] = thermodynamics.cid_to_dG0(cid)
        
        except MissingCompoundFormationEnergy:
            # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
            dG0_f[c] = pylab.nan
            ind_nan.append(c)
            html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td></tr>\n' % \
                              (kegg.cid2link(cid), cid, name))
    html_writer.write('</table>\n')
    
    bounds = [thermodynamics.bounds.get(cid, (None, None)) for cid in cids]
    res = {}
    try:
        res['pCr'] = find_pCr(S, dG0_f, c_mid=thermodynamics.c_mid, ratio=3.0, bounds=bounds)
        #res['PCR2'] = find_unfeasible_concentrations(S, dG0_f, c_range, c_mid=c_mid, bounds=bounds)
        res['MCMF'] = find_mcmf(S, dG0_f, c_range=thermodynamics.c_range, bounds=bounds)
    except LinProgNoSolutionException:
        html_writer.write('<b>No feasible solution found, cannot calculate the Margin</b>')
    
    # plot the profile graph
    pylab.rcParams['text.usetex'] = False
    pylab.rcParams['legend.fontsize'] = 10
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 12
    pylab.rcParams['lines.linewidth'] = 2
    pylab.rcParams['lines.markersize'] = 5
    pylab.rcParams['figure.figsize'] = [8.0, 6.0]
    pylab.rcParams['figure.dpi'] = 100

    # plot the thermodynamic profile in standard conditions
    
    profile_fig = pylab.figure()
    profile_fig.hold(True)

    pylab.title('Thermodynamic profile', figure=profile_fig)
    pylab.ylabel('cumulative dG [kJ/mol]', figure=profile_fig)
    pylab.xlabel('Reaction KEGG ID', figure=profile_fig)
    pylab.xticks(pylab.arange(1, Nr + 1), ['R%05d' % rids[i] for i in xrange(Nr)], fontproperties=FontProperties(size=8), rotation=30)

    dG0_r = pylab.zeros((Nr, 1))
    for r in range(Nr):
        reactants = pylab.find(S[r,:])
        dG0_r[r, 0] = pylab.dot(S[r, reactants], dG0_f[reactants])

    nan_indices = pylab.find(pylab.isnan(dG0_r))
    finite_indices = pylab.find(pylab.isfinite(dG0_r))
    if (len(nan_indices) > 0):
        dG0_r_finite = pylab.zeros((Nr, 1))
        dG0_r_finite[finite_indices] = dG0_r[finite_indices]
        cum_dG0_r = pylab.cumsum([0] + [dG0_r_finite[r, 0] * fluxes[r] for r in range(Nr)])
    else:
        cum_dG0_r = pylab.cumsum([0] + [dG0_r[r, 0] * fluxes[r] for r in range(Nr)])
    pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG0_r, figure=profile_fig, label='Standard [1M]')
    
    # plot the thermodynamic profile for the different optimization schemes
    
    pylab.grid(True, figure=profile_fig)
    for optimization in res.keys():
        (dG_f, conc, score) = res[optimization]
        if (score == None):
            continue

        dG_r = pylab.dot(S, dG_f)
        cum_dG_r = pylab.cumsum([0] + [dG_r[i, 0] * fluxes[i] for i in range(Nr)])
        pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG_r, figure=profile_fig, label='%s = %.1f' % (optimization, score))

    pylab.legend()
    html_writer.embed_matplotlib_figure(profile_fig, width=480, height=360)

    # plot the optimal metabolite concentrations for the different optimization schemes

    for optimization in res.keys():
        (dG_f, conc, score) = res[optimization]
        if (score == None):
            continue

        dG_r = pylab.dot(S, dG_f)
        conc[ind_nan] = thermodynamics.c_mid # give all compounds with unknown dG0_f the middle concentration value

        conc_fig = pylab.figure()
        conc_fig.suptitle('Concentrations (%s = %.1f)' % (optimization, score))
        pylab.xscale('log', figure=conc_fig)
        pylab.ylabel('Compound KEGG ID', figure=conc_fig)
        pylab.xlabel('Concentration [M]', figure=conc_fig)
        pylab.yticks(range(Nc, 0, -1), ["C%05d" % cid for cid in cids], fontproperties=FontProperties(size=8))
        pylab.plot(conc, range(Nc, 0, -1), '*b', figure=conc_fig)

        x_min = conc.min() / 10
        x_max = conc.max() * 10
        y_min = 0
        y_max = Nc + 1
        
        for c in range(Nc):
            pylab.text(conc[c, 0] * 1.1, Nc - c, kegg.cid2name(cids[c]), \
                       figure=conc_fig, fontsize=6, rotation=0)
            (b_low, b_up) = bounds[c]
            if (b_low == None):
                b_low = x_min
            if (b_up == None):
                b_up = x_max
            pylab.plot([b_low, b_up], [Nc - c, Nc - c], '-k', linewidth=0.4)

        if (optimization == 'pCr'):
            c_range_opt = pC_to_range(score, c_mid=thermodynamics.c_mid, ratio=3.0)
            pylab.axvspan(c_range_opt[0], c_range_opt[1], facecolor='g', alpha=0.3, figure=conc_fig)
        else:
            pylab.axvspan(thermodynamics.c_range[0], thermodynamics.c_range[1], facecolor='r', alpha=0.3, figure=conc_fig)
        pylab.axis([x_min, x_max, y_min, y_max], figure=conc_fig)
        html_writer.embed_matplotlib_figure(conc_fig, width=420, height=360)

    # write all the results in tables as well

    for optimization in res.keys():
        (dG_f, conc, score) = res[optimization]
        html_writer.write('<p>Biochemical Compound Formation Energies (%s = %.1f)<br>\n' % (optimization, score))
        html_writer.write('<table border="1">\n')
        html_writer.write('  ' + '<td>%s</td>'*5 % ("KEGG CID", "Compound Name", "Concentration [M]", "dG'0_f [kJ/mol]", "dG'_f [kJ/mol]") + '\n')
        for c in range(Nc):
            cid = cids[c]
            name = kegg.cid2name(cid)

            if (pylab.isnan(dG0_f[c, 0])):
                html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n' % \
                                  (kegg.cid2link(cid), cid, name, "N/A", "N/A", "N/A"))
            else:
                html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>%.2g</td><td>%.2f</td><td>%.2f</td></tr>\n' % \
                                  (kegg.cid2link(cid), cid, name, conc[c, 0], dG0_f[c, 0], dG_f[c, 0]))
        html_writer.write('</table></p>\n')

        html_writer.write('<p>Biochemical Reaction Energies (%s = %.1f)<br>\n' % (optimization, score))
        html_writer.write('<table border="1">\n')
        html_writer.write('  ' + '<td>%s</td>'*3 % ("KEGG RID", "dG'0_r [kJ/mol]", "dG'_r [kJ/mol]") + '\n')
        dG_r = pylab.dot(S, dG_f)
        for r in range(Nr):
            rid = rids[r]
            if (pylab.isnan(dG0_r[r, 0])):
                html_writer.write('<tr><td><a href="%s" title="%s">R%05d</a></td><td>%s</td><td>%.2f</td></tr>\n' % \
                                  (kegg.rid2link(rid), kegg.rid2name(rid), rid, "N/A", dG_r[r, 0]))
            else:
                html_writer.write('<tr><td><a href="%s" title="%s">R%05d</a></td><td>%.2f</td><td>%.2f</td></tr>\n' % \
                                  (kegg.rid2link(rid), kegg.rid2name(rid), rid, dG0_r[r, 0], dG_r[r, 0]))
        html_writer.write('</table></p>\n')
        
    return res

def test_all_modules():
    from pygibbs.groups import GroupContribution
    gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_test")
    gc.init()
    c_range = (1e-6, 1e-2)
    c_mid = 1e-3
    pH = 8
    I = 0.1
    T = 300
    map_cid = {201:2, 454:8} # CIDs that should be mapped to other CIDs because they are unspecific (like NTP => ATP)
    
    cids_with_missing_dG_f = set()
    
    f = open("../res/feasibility.csv", "w")
    csv_output = csv.writer(f)
    csv_output.writerow(("MID", "module name", "pH", "I", "T", "pCr", "MCMF"))
    for mid in sorted(gc.kegg().mid2rid_map.keys()):
        module_name = gc.kegg().mid2name_map[mid]
        try:
            (S, rids, fluxes, cids) = gc.kegg().get_module(mid)
        except kegg.KeggMissingModuleException as e:
            sys.stderr.write("WARNING: " + str(e) + "\n")
            continue
        (Nr, Nc) = S.shape
        for pH in [5, 6, 7, 8, 9]:
            for I in [0.0, 0.1, 0.2, 0.3, 0.4]:
                dG0_f = pylab.zeros((Nc, 1))
                bounds = []
                for c in range(Nc):
                    cid = map_cid.get(cids[c], cids[c])
                    try:
                        pmap = gc.cid2pmap(cid)
                        dG0_f[c] = gc.pmap_to_dG0(pmap, pH, I, T)
                    except MissingCompoundFormationEnergy as e:
                        if (cid not in cids_with_missing_dG_f):
                            sys.stderr.write("Setting the dG0_f of C%05d to NaN because: %s\n"\
                                             % (cid, str(e)))
                            cids_with_missing_dG_f.add(cid)
                        dG0_f[c] = pylab.nan
            
                bounds = [gc.kegg().cid2bounds.get(cid, (None, None)) for cid in cids]

                try:
                    (dG_f, concentrations, pCr) = find_pCr(S, dG0_f, c_mid=c_mid, ratio=3.0, bounds=bounds)
                except LinProgNoSolutionException:
                    sys.stderr.write("M%05d: Pathway is theoretically infeasible\n" % mid)
                    pCr = None

                try:
                    (dG_f, concentrations, MCMF) = find_mcmf(S, dG0_f, c_range=c_range, bounds=bounds)
                except LinProgNoSolutionException:
                    sys.stderr.write("M%05d: Pathway is theoretically infeasible\n" % mid)
                    MCMF = None
                
                csv_output.writerow([mid, module_name, pH, I, T, pCr, MCMF])
                        
    f.close()

def test_single_modules(mids):
    from pygibbs.groups import GroupContribution
    html_writer = HtmlWriter("../res/thermodynamic_module_analysis.html")
    gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_test")
    gc.init()
    
    for mid in mids:
        html_writer.write("<h2>M%05d</h2>\n" % mid)
        (S, rids, fluxes, cids) = gc.kegg().get_module(mid)
        thermodynamic_pathway_analysis(S, rids, fluxes, cids, gc, gc.kegg(), html_writer)

if (__name__ == "__main__"):
    test_single_modules([5, 305, 719])
    #test_all_modules()