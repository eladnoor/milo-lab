import alberty, hatzimanikatis, sys
from nist import Nist
from pylab import *
from alberty import Alberty, MissingCompoundFormationEnergy
from hatzimanikatis import Hatzi
from common import *
from groups import GroupContribution, GroupDecompositionError
from nist_train import GradientAscent
from matplotlib import font_manager

A = Alberty()
H = Hatzi()
gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_test")
gc.init()
nist = Nist(gc.kegg())
grad = GradientAscent(gc)
grad.load_energies()

def pH_dependence():
    
    analyze_this_reaction = []
    I_mid = []
    I_tolerance = []
    T_mid = []
    T_tolerance = []
    
    analyze_this_reaction += [{2:-1, 31:-1, 8:1, 92:1}] # Glucose kinase
    I_mid += [0.01]
    I_tolerance += [0.02]
    T_mid += [303.1]
    T_tolerance += [0.1]
    
    analyze_this_reaction += [{1:-1, 1005:-1, 9:1, 65:1}] # L-serine kinase
    I_mid += [0.25]
    I_tolerance += [0.01]
    T_mid += [311.15]
    T_tolerance += [0.1]
    
    analyze_this_reaction += [{1:-1, 13:-1, 9:2}] # PPi hydrolase
    I_mid += [0.05]
    I_tolerance += [0.05]
    T_mid += [298.1]
    T_tolerance += [15]
    
    analyze_this_reaction += [{354:-1, 111:1, 118:1}] # Fructose bisphosphate aldolase
    I_mid += [0.01]
    I_tolerance += [0.011]
    T_mid += [300]
    T_tolerance += [12]
    
    figure()
    rcParams['text.usetex'] = True
    rcParams['legend.fontsize'] = 4
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.size'] = 6
    rcParams['lines.linewidth'] = 0.5
    rcParams['lines.markersize'] = 3       
    for i in range(len(analyze_this_reaction)):
        subplot(2,2,i+1)
        
        sys.stdout.write("Compound parameters from Alberty's table:\n")
        for cid in analyze_this_reaction[i].keys():
            A.display_pmap(cid)
        
        sys.stdout.write("Compound parameters from Hatzimanikatis' table:\n")
        for cid in analyze_this_reaction[i].keys():
            H.display_pmap(cid)
        
        M_obs = []
        sys.stdout.write("%5s | %5s | %5s | %6s | %6s | %s\n" % ("Match", "pH", "I", "T", "dG0(N)", "link"))
        for row in nist.data:
            if (analyze_this_reaction[i] != None and analyze_this_reaction[i] != row[6]):
                continue
            try:
                sparse_reaction = row[6]
                #evaluation = row[3] # A, B, C, D
                [Keq, T, I, pH] = row[8:12]
                dG0_N = -R*T*log(Keq)
                if (analyze_this_reaction[i] == sparse_reaction and abs(I-I_mid[i])<I_tolerance[i] and abs(T-T_mid[i])<T_tolerance[i]):
                    M_obs.append([pH, dG0_N])
                    sys.stdout.write(" ***  | ")
                else:
                    sys.stdout.write("      | ")
                sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f | %s\n" % (pH, I, T, dG0_N, row[0]))
            except MissingCompoundFormationEnergy:
                continue
        if (len(M_obs) == 0):
            sys.stderr.write("There are now data points matching this reaction with the specifice I and T")
            continue
        M_obs = matrix(M_obs)
        
        leg = ['NIST']
        
        pH_range = arange(M_obs[:,0].min()-1, M_obs[:,0].max()+1, 0.01)
        M_est = []
        I_low  = max(I_mid[i]-I_tolerance[i], 0.0)
        I_high = I_mid[i]+I_tolerance[i]
        for pH in pH_range:
            predictions = []
            for predictor in [A, H, grad]:
                predictions.append(predictor.reaction_to_dG0(sparse_reaction, pH, I=I_low, T=T_mid[i]))
                predictions.append(predictor.reaction_to_dG0(sparse_reaction, pH, I=I_mid[i], T=T_mid[i]))
                predictions.append(predictor.reaction_to_dG0(sparse_reaction, pH, I=I_high, T=T_mid[i]))
            M_est.append(predictions)
        M_est = matrix(M_est)
        
        plot(M_obs[:,0], M_obs[:,1], 'co')
        
        plot(pH_range, M_est[:,0], 'b:')
        plot(pH_range, M_est[:,1], 'b-')
        plot(pH_range, M_est[:,2], 'b:')
    
        plot(pH_range, M_est[:,3], 'g:')
        plot(pH_range, M_est[:,4], 'g-')
        plot(pH_range, M_est[:,5], 'g:')
    
        plot(pH_range, M_est[:,6], 'r:')
        plot(pH_range, M_est[:,7], 'r-')
        plot(pH_range, M_est[:,8], 'r:')
    
        for predictor_name in ['Alberty', 'Hatzi', 'Rugged']:
            for I in [I_low, I_mid[i], I_high]: 
                leg.append("%s (I = %.2f)" % (predictor_name, I))
        legend(leg, loc='best')
        
        if (i == 2 or i == 3):
            xlabel('pH')
        if (i == 0 or i == 2):
            ylabel(r"$\Delta_r G'^\circ$ [kJ/mol]")
        s_title = gc.kegg().sparse_reaction_to_string(sparse_reaction, cids=False) + "\n"
        s_title += "($I = %.2f \pm %.2f$ $M$, $T = %.1f \pm %.1f$ $K$)" % (I_mid[i], I_tolerance[i], T_mid[i]-273.15, T_tolerance[i])
        title(s_title, fontsize=5)
    savefig('../res/compare_pH.pdf', format='pdf')

def map_rid_to_nist_rowids():
    rid_to_nist_rowids = {}
    for rid in gc.kegg().get_all_rids():
        sparse_reaction = gc.kegg().rid2sparse_reaction(rid)
        for rowid in range(len(nist.data)):
            row = nist.data[rowid]
            if (sparse_reaction == row[6]):
                rid_to_nist_rowids[rid] = rid_to_nist_rowids.get(rid, []) + [rowid]
    return rid_to_nist_rowids

def calculate_uncertainty(sparse_reaction, min_C, max_C, T):
    N_subs = 0
    N_prod = 0
    for (cid, coeff) in sparse_reaction.iteritems():
        if (cid == 1 or cid == 80): # uncertainty in H2O and H+ are ignored
            continue
        if (coeff > 0):
            N_prod += coeff
        else:
            N_subs -= coeff
    ddG_min = R * T * (N_prod * log(min_C) - N_subs * log(max_C))
    ddG_max = R * T * (N_prod * log(max_C) - N_subs * log(min_C))
    return (ddG_min, ddG_max)

def uncertainty_comparison():
    """
        Plots the uncertainty in the dG of each reaction in KEGG, against the RMSE of the predictions (compared to the NIST measurements).
        The x-value is a function of the number of substrates and number of products of each reaction (excluding H2O and H+).
        The y-value of each dot is calculated by the RMSE of the observation vs. estimation across all measurements of the same reaction.
    """
    limits = [(-5, -3), (-5, -2), (-6, -2)] # in log10 scale
    rcParams['text.usetex'] = True
    rcParams['legend.fontsize'] = 8
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.size'] = 8
    rcParams['lines.linewidth'] = 0.4
    rcParams['lines.markersize'] = 3

    figure(figsize=(12,3.5))
    for i in range(len(limits)):
        (min_C, max_C) = limits[i]
        
        rid_to_nist_rowids = map_rid_to_nist_rowids()
        data_mat = []
        for (rid, rowids) in rid_to_nist_rowids.iteritems():
            sparse_reaction = gc.kegg().rid2sparse_reaction(rid)
            try:
                error_mat = []
                for rowid in rowids:
                    row = nist.data[rowid]
                    #evaluation = row[3] # A, B, C, D
                    [Keq, T, I, pH] = row[8:12]
                    dG0_obs = -R*T*log(Keq)
                    dG0_est = [predictor.reaction_to_dG0(sparse_reaction, pH, I, T) for predictor in [A, H, grad]]
                    error_mat.append([(dG0_obs - x) for x in dG0_est])
                error_mat = array(error_mat)
                rmse = sqrt(mean(error_mat**2, 0))
                (ddG_min, ddG_max) = calculate_uncertainty(sparse_reaction, min_C=10**min_C, max_C=10**max_C, T=300)
                data_mat.append([ddG_max - ddG_min, rmse[0], rmse[1], rmse[2]])
            except MissingCompoundFormationEnergy:
                continue
        
        data_mat = matrix(data_mat)
    
        subplot(1,len(limits),i+1)
        hold(True)
        plot(data_mat[:,0], data_mat[:,1:], '.')
        plot([0, 200], [0, 200], '--k')
        axis('scaled')
        xlabel(r"uncertainty in $\Delta_r G$ due to concentrations [kJ/mol]")
        if (i == 0):
            ylabel(r"RMSE of $\Delta_r G^\circ$ estimation [kJ/mol]")
        title(r"$10^{%g}M$ $<$ [c] $<$ $10^{%g}M$" % (min_C, max_C))
        legend(['Alberty', 'Hatzimanikatis', 'Rugged'], loc="upper left")

    savefig('../res/compare_uncertainty.pdf', format='pdf')
    
if (__name__ == "__main__"):
    #pH_dependence()
    uncertainty_comparison()