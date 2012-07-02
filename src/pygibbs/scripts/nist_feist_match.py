from pygibbs.nist import Nist
from pygibbs.feist_ecoli import Feist
import csv
import numpy as np
from pygibbs.dissociation_constants import DissociationConstants,\
    MissingDissociationConstantError


(pH, I, pMg, T) = (7.0, 0.1, 14, 298.15)

nist_csv = csv.writer(open('../res/nist_rt_data.tsv', 'w'), delimiter='\t')
nist_csv.writerow(('URL', "dG'0", 'pH', 'I', 'T', 'dG0', 'ddG0', 'match in iAF1260'))

feist_csv = csv.writer(open('../res/nist_feist_match.tsv', 'w'), delimiter='\t')
feist_csv.writerow(('name', "dG'0", 'pH', 'I', 'T', 'dG0', 'ddG0'))

nist = Nist()
dissociation = DissociationConstants.FromPublicDB()
cid2nH_nMg = dissociation.GetCid2nH_nMg(pH, I, pMg, T)

for i, r in enumerate(Feist.FromFiles().reactions):
    if r.name in ['NTPP9', 'DHQS', 'AICART']: # these reactions were not in NIST when the UGCM model was created
        continue
    
    nist_rows = nist.SelectRowsFromNist(r, check_reverse=True)
    if nist_rows == []:
        continue

    try:
        print r.name
        
        dG0_list = []
        for row in nist_rows:
            if row.reaction == r:
                dG0_prime = row.dG0_r
            elif row.reaction.reverse() == r:
                dG0_prime = -row.dG0_r
            else:
                raise Exception("internal error")
            
            ddG0 = dissociation.ReverseTransformReaction(r, 
                      row.pH, row.I, row.pMg, row.T, cid2nH_nMg=cid2nH_nMg,
                      suppress_missing_pka_exception=False)
            dG0 = dG0_prime - ddG0
            nist_csv.writerow((row.url, "%.1f" % dG0_prime, row.pH, row.I, row.T, "%.1f" % dG0, "%.1f" % ddG0, r.name))
            dG0_list.append(dG0)
        
        ddG0_avg = dissociation.ReverseTransformReaction(r,
                      pH, I, pMg, T, cid2nH_nMg=cid2nH_nMg,
                      suppress_missing_pka_exception=False)
        dG0_avg = np.mean(dG0_list)
        dG0_prime = dG0_avg + ddG0_avg
        feist_csv.writerow((r.name, "%.1f" % dG0_prime, pH, I, T, "%.1f" % dG0_avg, "%.1f" % ddG0_avg))
    except MissingDissociationConstantError:
        continue