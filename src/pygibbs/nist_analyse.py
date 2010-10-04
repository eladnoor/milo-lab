import csv
from nist import Nist
import kegg

KEGG = kegg.Kegg()

def cid2name(cid, KEGG):
    return "\"" + KEGG.cid2name(cid) + "\""

def load_cid_set(train_csv_fname):
    """
        Read the training data from a CSV file
    """
    csv_reader = csv.reader(open(train_csv_fname))
    csv_reader.next()
    cid_set = set()
    for row in csv_reader:
        #(smiles, cid, compoud_name, dG0, dH0, z, nH, Mg, use_for, ref, remark) = row
        if (row[8] in ['skip']):
            continue
        cid = int(row[1])
        if (cid > 0):
            cid_set.add(cid)
    return cid_set
        
#############################################################################################################

nist = Nist(KEGG)
known_cids = load_cid_set('../data/dG0_seed.csv')

pairs = set()
dot_file = open('../res/nist.dot', 'w')
dot_file.write('graph G {\n')

one_step_cids = set()
coupled_cids = set()

for row in nist.data:
    sparse_reaction = row[6]
    unknown_cids = list(set(sparse_reaction.keys()).difference(known_cids))
    if (len(unknown_cids) == 1):
        one_step_cids.add(unknown_cids[0])
    elif (len(unknown_cids) == 2):
        coupled_cids.add((min(unknown_cids), max(unknown_cids)))

for cid in one_step_cids:
    dot_file.write("%s [fillcolor = red]\n" % cid2name(cid, KEGG))

for (cid1, cid2) in coupled_cids:
    dot_file.write("%s -- %s\n" % (cid2name(cid1, KEGG), cid2name(cid2, KEGG)))

dot_file.write('}\n')
