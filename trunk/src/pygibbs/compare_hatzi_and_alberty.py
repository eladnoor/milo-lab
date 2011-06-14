from nist import Nist
import pylab
from alberty import Alberty, MissingCompoundFormationEnergy
from hatzimanikatis import Hatzi
from groups import GroupContribution
from thermodynamic_constants import R
import logging
import sys
from toolbox.database import SqliteDatabase
from pygibbs.kegg_reaction import Reaction

A = Alberty()
H = Hatzi()
db = SqliteDatabase('../res/gibbs.sqlite')
gc = GroupContribution(db)
gc.init()
nist = Nist()

def pH_dependence():
    
    analyze_this_reaction = []
    I_mid = []
    I_tolerance = []
    T_mid = []
    T_tolerance = []
    
    analyze_this_reaction += [Reaction(['glucose kinase'], {2:-1, 31:-1, 8:1, 92:1})]
    I_mid += [0.01]
    I_tolerance += [0.02]
    T_mid += [303.1]
    T_tolerance += [0.1]
    
    analyze_this_reaction += [Reaction(['L-serine kinase'], {1:-1, 1005:-1, 9:1, 65:1})]
    I_mid += [0.25]
    I_tolerance += [0.01]
    T_mid += [311.15]
    T_tolerance += [0.1]
    
    analyze_this_reaction += [Reaction(['pyrophosphatase'], {1:-1, 13:-1, 9:2}, rid=4)]
    I_mid += [0.05]
    I_tolerance += [0.05]
    T_mid += [298.1]
    T_tolerance += [15]
    
    analyze_this_reaction += [Reaction(['Fructose bisphosphate aldolase'], {354:-1, 111:1, 118:1})]
    I_mid += [0.01]
    I_tolerance += [0.011]
    T_mid += [300]
    T_tolerance += [12]
    
    pylab.figure()
    pylab.rcParams['text.usetex'] = True
    pylab.rcParams['legend.fontsize'] = 4
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 6
    pylab.rcParams['lines.linewidth'] = 0.5
    pylab.rcParams['lines.markersize'] = 3       
    for i, reaction in enumerate(analyze_this_reaction):
        pylab.subplot(2,2,i+1)
        
        logging.info("Compound parameters from Alberty's table:")
        for cid in reaction.get_cids():
            A.cid2PseudoisomerMap(cid).Display()
        
        logging.info("Compound parameters from Hatzimanikatis' table:")
        for cid in reaction.get_cids():
            H.cid2PseudoisomerMap(cid).Display()
        
        M_obs = []
        sys.stdout.write("%5s | %5s | %5s | %6s | %6s | %s\n" % ("Match", "pH", "I", "T", "dG0(N)", "link"))
        for row in nist.data:
            if (reaction != None and reaction != row.reaction):
                continue
            try:
                #evaluation = row[3] # A, B, C, D
                if (reaction == row.reaction and 
                    abs(row.I-I_mid[i]) < I_tolerance[i] and 
                    abs(row.T-T_mid[i]) < T_tolerance[i]):
                    M_obs.append([row.pH, row.dG0_r])
                    sys.stdout.write(" ***  | ")
                else:
                    sys.stdout.write("      | ")
                sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f | %s\n" % (row.pH, row.I, row.T, row.dG0_r, row.ref_id))
            except MissingCompoundFormationEnergy:
                continue
        if len(M_obs) == 0:
            sys.stderr.write("There are now data points matching this reaction with the specific I and T")
            continue
        M_obs = pylab.matrix(M_obs)
        
        leg = ['NIST']
        
        pH_range = pylab.arange(M_obs[:,0].min()-1, M_obs[:,0].max()+1, 0.01)
        M_est = []
        I_low  = max(I_mid[i]-I_tolerance[i], 0.0)
        I_high = I_mid[i]+I_tolerance[i]
        for pH in pH_range:
            predictions = []
            for predictor in [A, H]:
                predictions.append(reaction.PredictReactionEnergy(predictor, pH=pH, I=I_low ,T=T_mid[i]))
                predictions.append(reaction.PredictReactionEnergy(predictor, pH=pH, I=I_mid[i] ,T=T_mid[i]))
                predictions.append(reaction.PredictReactionEnergy(predictor, pH=pH, I=I_high ,T=T_mid[i]))
            M_est.append(predictions)
        M_est = pylab.matrix(M_est)
        
        pylab.plot(M_obs[:,0], M_obs[:,1], 'co')
        
        pylab.plot(pH_range, M_est[:,0], 'b:')
        pylab.plot(pH_range, M_est[:,1], 'b-')
        pylab.plot(pH_range, M_est[:,2], 'b:')
    
        pylab.plot(pH_range, M_est[:,3], 'g:')
        pylab.plot(pH_range, M_est[:,4], 'g-')
        pylab.plot(pH_range, M_est[:,5], 'g:')
    
        pylab.plot(pH_range, M_est[:,6], 'r:')
        pylab.plot(pH_range, M_est[:,7], 'r-')
        pylab.plot(pH_range, M_est[:,8], 'r:')
    
        for predictor_name in ['Alberty', 'Hatzi', 'Rugged']:
            for I in [I_low, I_mid[i], I_high]: 
                leg.append("%s (I = %.2f)" % (predictor_name, I))
        pylab.legend(leg, loc='best')
        
        if (i == 2 or i == 3):
            pylab.xlabel('pH')
        if (i == 0 or i == 2):
            pylab.ylabel(r"$\Delta_r G'^\circ$ [kJ/mol]")
        s_title = gc.kegg.reaction2string(reaction, cids=False) + "\n"
        s_title += "($I = %.2f \pm %.2f$ $M$, $T = %.1f \pm %.1f$ $K$)" % (I_mid[i], I_tolerance[i], T_mid[i]-273.15, T_tolerance[i])
        pylab.title(s_title, fontsize=5)
    pylab.savefig('../res/compare_pH.pdf', format='pdf')

def map_rid_to_nist_rowids():
    rid_to_nist_rowids = {}
    for rid in gc.kegg.get_all_rids():
        reaction = gc.kegg.rid2reaction(rid)
        for rowid in range(len(nist.data)):
            row = nist.data[rowid]
            if reaction == row.reaction:
                rid_to_nist_rowids[rid] = rid_to_nist_rowids.get(rid, []) + [rowid]
    return rid_to_nist_rowids

def calculate_uncertainty(reaction, min_C, max_C, T):
    N_subs = 0
    N_prod = 0
    for cid, coeff in reaction.sparse.iteritems():
        if cid == 1 or cid == 80: # uncertainty in H2O and H+ are ignored
            continue
        if coeff > 0:
            N_prod += coeff
        else:
            N_subs -= coeff
    ddG_min = R * T * (N_prod * pylab.log(min_C) - N_subs * pylab.log(max_C))
    ddG_max = R * T * (N_prod * pylab.log(max_C) - N_subs * pylab.log(min_C))
    return (ddG_min, ddG_max)

def uncertainty_comparison():
    """
        Plots the uncertainty in the dG of each reaction in KEGG, against the RMSE of the predictions (compared to the NIST measurements).
        The x-value is a function of the number of substrates and number of products of each reaction (excluding H2O and H+).
        The y-value of each dot is calculated by the RMSE of the observation vs. estimation across all measurements of the same reaction.
    """
    limits = [(-5, -3), (-5, -2), (-6, -2)] # in log10 scale
    pylab.rcParams['text.usetex'] = True
    pylab.rcParams['legend.fontsize'] = 8
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 8
    pylab.rcParams['lines.linewidth'] = 0.4
    pylab.rcParams['lines.markersize'] = 3

    pylab.figure(figsize=(12,3.5))
    for i in range(len(limits)):
        (min_C, max_C) = limits[i]
        
        rid_to_nist_rowids = map_rid_to_nist_rowids()
        data_mat = []
        for (rid, rowids) in rid_to_nist_rowids.iteritems():
            reaction = gc.kegg.rid2reaction(rid)
            try:
                error_mat = []
                for rowid in rowids:
                    row = nist.data[rowid]
                    #evaluation = row[3] # A, B, C, D
                    dG0_est = [reaction.PredictReactionEnergy(predictor, pH=row.pH, I=row.I, T=row.T)
                               for predictor in [A, H]]
                    error_mat.append([(row.dG0_r - x) for x in dG0_est])
                error_mat = pylab.array(error_mat)
                rmse = pylab.sqrt(pylab.mean(error_mat**2, 0))
                (ddG_min, ddG_max) = calculate_uncertainty(reaction, min_C=10**min_C, max_C=10**max_C, T=300)
                data_mat.append([ddG_max - ddG_min, rmse[0], rmse[1], rmse[2]])
            except MissingCompoundFormationEnergy:
                continue
        
        data_mat = pylab.matrix(data_mat)
    
        pylab.subplot(1,len(limits),i+1)
        pylab.hold(True)
        pylab.plot(data_mat[:,0], data_mat[:,1:], '.')
        pylab.plot([0, 200], [0, 200], '--k')
        pylab.axis('scaled')
        pylab.xlabel(r"uncertainty in $\Delta_r G$ due to concentrations [kJ/mol]")
        if (i == 0):
            pylab.ylabel(r"RMSE of $\Delta_r G^\circ$ estimation [kJ/mol]")
        pylab.title(r"$10^{%g}M$ $<$ [c] $<$ $10^{%g}M$" % (min_C, max_C))
        pylab.legend(['Alberty', 'Hatzimanikatis', 'Rugged'], loc="upper left")

    pylab.savefig('../res/compare_uncertainty.pdf', format='pdf')
    
if (__name__ == "__main__"):
    #pH_dependence()
    uncertainty_comparison()