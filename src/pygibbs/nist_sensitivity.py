from groups import GroupContribution
from nist_train import GradientAscent
from alberty import Alberty
from nist import Nist
from copy import deepcopy
import csv, sys

################################################################################

def sensitivity_analysis_for_gradient_ascent(gc, nist, cid2pmap, max_i=250, n_begin=0):
    # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
    # Use DETERMINISTIC gradient ascent
    # run the simulation many times, each time removing a fraction of the reactions from the training set, and
    # evaluating the MSE only on that fraction (i.e. it will be the test set)
    
    grad = GradientAscent(gc)
    grad.cid2pmap = deepcopy(cid2pmap)
    grad.load_nist_data(nist, skip_missing_reactions=True)
    
    res_file = open('../res/leave1out_report.csv', 'w')
    csv_results = csv.writer(res_file)
    csv_results.writerow(["N", "dG0_obs", "dG0_est_before", "dG0_est_after", "iterations", "reaction", "pH", "I", "T", "evaluation", "no. measurements"])
    
    N = len(grad.data)
    for n in range(n_begin, N):
        (sparse_reaction, pH, I, T, evaluation, dG0_obs) = grad.data[n]
        n_measurements = min([nist.cid2count[cid] for cid in sparse_reaction.keys()])
        reaction_str = gc.kegg().sparse_reaction_to_string(sparse_reaction, cids=True)
        sys.stderr.write("%04d) %s, %.1f, %.2f, %.1f, %s, %d, " % (n, reaction_str, pH, I, T, evaluation, n_measurements))

        dG0_est_before = grad.reaction_to_dG0(sparse_reaction, pH, I, T)
        sys.stderr.write("E-before = %.2f, " % abs(dG0_obs - dG0_est_before))

        grad.cid2pmap = deepcopy(cid2pmap)
        grad.train_rowids = range(n) + range(n+1, N)
        grad.update_cache()
        grad.deterministic_hill_climb(max_i=max_i, verbose=False)
        dG0_est_after = grad.reaction_to_dG0(sparse_reaction, pH, I, T)
        
        csv_results.writerow([n, dG0_obs, dG0_est_before, dG0_est_after, reaction_str, max_i, pH, I, T, evaluation, n_measurements])
        res_file.flush()
        sys.stderr.write("E-after = %.2f\n" % abs(dG0_obs - dG0_est_after))
    
    #plot(dG0_obs_vec, dG0_est_after_vec, '.')
    #r2 = GradientAscent.calc_r2(dG0_obs_vec, dG0_est_after_vec)
    #title(r'N = %d, r$^2$ = %.2f' % (len(dG0_obs_vec), r2), fontsize=14)
    #xlabel(r'$\Delta_{obs} G^\circ$ [kJ/mol]', fontsize=14)
    #ylabel(r'$\Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
    #min_x = min(dG0_obs_vec)
    #max_x = max(dG0_obs_vec)
    #plot([min_x, max_x], [min_x, max_x], 'k--')
    #axis([-60, 60, -60, 60])
    #savefig("leave1out_meas_vs_obs.pdf", format='pdf')

def evaluate(gc, nist, cid2pmap):
    # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
    # Use DETERMINISTIC gradient ascent
    # run the simulation many times, each time removing a fraction of the reactions from the training set, and
    # evaluating the MSE only on that fraction (i.e. it will be the test set)
    
    grad = GradientAscent(gc)
    grad.cid2pmap = deepcopy(cid2pmap)
    grad.load_nist_data(nist, skip_missing_reactions=True)

    res_file = open('../res/evaluation_report.csv', 'w')
    csv_results = csv.writer(res_file)
    csv_results.writerow(["N", "dG0_obs", "dG0_est", "reaction", "pH", "I", "T", "evaluation"])
    
    N = len(grad.data)
    for n in range(n_begin, N):
        (sparse_reaction, pH, I, T, evaluation, dG0_obs) = grad.data[n]
        n_measurements = min([nist.cid2count[cid] for cid in sparse_reaction.keys()])
        reaction_str = gc.kegg().sparse_reaction_to_string(sparse_reaction, cids=True)
        dG0_est = grad.reaction_to_dG0(sparse_reaction, pH, I, T)
        csv_results.writerow([n, dG0_obs, dG0_est, reaction_str, pH, I, T, evaluation, n_measurements])
        res_file.flush()
            
################################################################################

if (len(sys.argv) > 1):
    n_begin = int(sys.argv[1])
else:
    n_begin = 0

gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_test")
gc.init()
nist = Nist(gc.kegg())
alberty = Alberty()
sensitivity_analysis_for_gradient_ascent(gc, nist, alberty.cid2pmap, max_i=250, n_begin=n_begin)
#evaluate(gc, nist, alberty.cid2pmap)