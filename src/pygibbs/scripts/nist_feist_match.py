from pygibbs.nist import Nist
from pygibbs.feist_ecoli import Feist
import csv

csv_writer = csv.writer(open('../res/nist_feist_match.tsv', 'w'), delimiter='\t')
csv_writer.writerow(('name', "dG'0", 'pH', 'I', 'T', 'evaluation', 'reference ID'))
nist = Nist()
for i, r in enumerate(Feist.FromFiles().reactions):
    print r.name
    for row in nist.SelectRowsFromNist(r, check_reverse=True):
        if row.reaction == r:
            dG0_prime = row.dG0_r
        else:
            dG0_prime = -row.dG0_r
        csv_writer.writerow((r.name, dG0_prime, row.pH, row.I, row.T, row.evaluation, row.ref_id))
    