import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rms_flat
from pygibbs.feist_ecoli import Feist
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.nist import Nist
from pygibbs.unified_group_contribution import UnifiedGroupContribution
from toolbox.database import SqliteDatabase
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg import Kegg
from pygibbs.thermodynamic_constants import R
import logging
import csv

DATA_FNAME = '../res/compare_ugcm_to_hatzi.txt'
CSV_FNAME = '../res/compare_ugcm_to_hatzi.csv'
FIG_FNAME = '../res/compare_ugcm_to_iAF1260'
pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)
pH_min, pH_max = (6.0, 8.0)

if True:
    estimators = LoadAllEstimators()
    nist = Nist()
    feist = Feist.FromFiles()
    
    reactions = {}
    
    csv_writer = csv.DictWriter(open(CSV_FNAME, 'w'),
        fieldnames = ['name', 'reaction', 'E[dG''0] nist', 'dG''0 iAF1260', 'dG''0 UGCM'],
        extrasaction='ignore')
    csv_writer.writeheader()

    data = []
    for i in xrange(len(feist.reactions)):
        r = feist.reactions[i]
        logging.info("Reaction %s: %s" % (r.name, r.FullReactionString()))
        dG0_feist = feist.dG0s[i]
        
        nist_rows = nist.SelectRowsFromNist(r, check_reverse=True)
        dG0_list = []
        for row_data in nist_rows:
            if pH_min < row_data.pH < pH_max:
                if row_data.reaction == r:
                    dG0_list.append(row_data.dG0_r)
                else:
                    dG0_list.append(-row_data.dG0_r)
        if len(dG0_list) > 0:
            dG0_nist = np.mean(dG0_list)
        else:
            dG0_nist = np.nan
            
        d = {}
        d['name'] = r.name
        d['reaction'] = r.FullReactionString()
        d['E[dG''0] nist'] = dG0_nist
        dG0_prime_hatzi = np.nan
        dG0_prime_ugcm = np.nan
        try:
            dG0_prime_hatzi = r.PredictReactionEnergy(estimators['hatzi_gc'], pH=pH, pMg=pMg, I=I, T=T)
            dG0_prime_ugcm = r.PredictReactionEnergy(estimators['UGC'], pH=pH, pMg=pMg, I=I, T=T)
        except ValueError:
            pass

        d['dG''0 iAF1260'] = dG0_prime_hatzi
        d['dG''0 UGCM'] = dG0_prime_ugcm
        data.append((dG0_prime_hatzi, dG0_feist, dG0_prime_ugcm, dG0_nist))
            
        csv_writer.writerow(d)

    np.savetxt(DATA_FNAME, np.matrix(data), fmt='%.2f', delimiter=',', newline='\n')

data = np.loadtxt(DATA_FNAME, dtype='float', delimiter=',')

#plt.plot(data[:, 0], data[:, 1], '.')

feist_idx = set(np.nonzero(np.isfinite(data[:, 1]))[0].flat)
ugcm_idx = set(np.nonzero(np.isfinite(data[:, 2]))[0].flat)
nist_idx = set(np.nonzero(np.isfinite(data[:, 3]))[0].flat)

comp_idx = list(feist_idx.intersection(ugcm_idx).intersection(nist_idx))
minG, maxG = (np.min(data[comp_idx, 0]), np.max(data[comp_idx, 0]))

plt.figure(figsize=(10, 5), dpi=90)
plt.subplot(1,2,1)
err_feist_nist = data[comp_idx, 1] - data[comp_idx, 3]
rms_feist_nist = rms_flat(err_feist_nist)
plt.plot(data[comp_idx, 1], data[comp_idx, 3], '.g')
plt.plot([minG, maxG], [minG, maxG], ':k')
plt.ylabel('TECRDB observation [kJ/mol]')
plt.xlabel('value in iAF1260 [kJ/mol]')
plt.title('N = %d, RMSE = %.1f [kJ/mol]' % (len(comp_idx), rms_feist_nist))

plt.subplot(1,2,2)
err_ugcm_nist = data[comp_idx, 2] - data[comp_idx, 3]
rms_ugcm_nist = rms_flat(err_ugcm_nist)
plt.plot(data[comp_idx, 2], data[comp_idx, 3], '.g')
plt.plot([minG, maxG], [minG, maxG], ':k')
plt.ylabel('TECRDB observation [kJ/mol]')
plt.xlabel('UGCM estimation [kJ/mol]')
plt.title('N = %d, RMSE = %.1f [kJ/mol]' % (len(comp_idx), rms_ugcm_nist))
plt.tight_layout()
plt.savefig(FIG_FNAME + "_fig1.svg", fmt='.svg')

plt.figure(figsize=(10, 5), dpi=90)
plt.subplot(1,2,1)
err_feist_ugcm1 = data[comp_idx, 1] - data[comp_idx, 2]
rms_feist_ugcm1 = rms_flat(err_feist_ugcm1)
plt.plot(data[comp_idx, 1], data[comp_idx, 2], '.g')
plt.plot([minG, maxG], [minG, maxG], ':k')
plt.xlabel('value in iAF1260 [kJ/mol]')
plt.ylabel('UGCM estimation [kJ/mol]')
plt.title('Observed data, N = %d, RMSE = %.1f [kJ/mol]' % (len(comp_idx), rms_feist_ugcm1))

plt.subplot(1,2,2)
non_nist_idx = list(feist_idx.intersection(ugcm_idx).difference(nist_idx))
err_feist_ugcm2 = data[non_nist_idx, 1] - data[non_nist_idx, 2]
rms_feist_ugcm2 = rms_flat(err_feist_ugcm2)
plt.plot(data[non_nist_idx, 1], data[non_nist_idx, 2], '.g')
plt.plot([minG, maxG], [minG, maxG], ':k')
plt.xlabel('value in iAF1260 [kJ/mol]')
plt.ylabel('UGCM estimation [kJ/mol]')
plt.title('Unobserved data, N = %d, RMSE = %.1f [kJ/mol]' % (len(non_nist_idx), rms_feist_ugcm2))
plt.tight_layout()
plt.savefig(FIG_FNAME + "_fig2.svg", fmt='.svg')

plt.figure(figsize=(6, 6), dpi=90)
bins = np.arange(-30, 30, 2)
plt.hist([err_feist_nist, err_ugcm_nist], bins=bins, histtype='bar', cumulative=False, normed=False)
plt.xlabel('Error in kJ/mol')
plt.ylabel('# of reactions')
plt.legend(['value in iAF1260', 'UGCM estimation'])
plt.savefig(FIG_FNAME + "_fig3.svg", fmt='.svg')

db = SqliteDatabase('../res/gibbs.sqlite', 'w')
ugc = UnifiedGroupContribution(db)
ugc.LoadGroups(True)
ugc.LoadObservations(True)
ugc.LoadGroupVectors(True)
ugc.LoadData(True)
ugc.init()
r_list = []
#r_list += [Reaction.FromFormula("C00036 + C00044 = C00011 + C00035 + C00074")]
#r_list += [Reaction.FromFormula("C00003 + C00037 + C00101 = C00004 + C00011 + C00014 + C00080 + C00143")] # glycine synthase
r_list += [Reaction.FromFormula("C00001 + C00002 + C00064 + C04376 => C00008 + C00009 + C00025 + C04640")]
#r_list += [Reaction.FromFormula("C00001 + 2 C00002 + C00064 + C00288 <=> 2 C00008 + C00009 + C00025 + C00169")]


kegg = Kegg.getInstance()
S, cids = kegg.reaction_list_to_S(r_list)

logging.getLogger('').setLevel(logging.DEBUG)

dG0_prime = ugc.GetTransfromedReactionEnergies(S, cids, pH=pH, I=I, pMg=pMg, T=T)
RT = R * T
for i in xrange(len(r_list)):
    r_list[i].Balance()
    dG0_r, parts, dG0_r_pgc = ugc.GetChemicalReactionEnergies(S[:, i], cids)
    print r_list[i].FullReactionString(show_cids=False)
    print ('UGC: dG0 = %.1f = ' % dG0_r.sum(0)) + ' + '.join('%.1f' % d for d in dG0_r.flat)
    print ('PGC: dG0 = %.1f = ' % dG0_r_pgc.sum(0)) + ' + '.join('%.1f' % d for d in dG0_r_pgc.flat)
    print "dG0' = %.1f" % (dG0_prime[0, i])
