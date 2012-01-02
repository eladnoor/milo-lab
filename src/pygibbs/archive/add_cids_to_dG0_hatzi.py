import kegg
import csv
import sys

###########################################################################################
def find_cid_in_kegg(kegg, name_list):
    for name in name_list:
        if (kegg.name2cid(name) != None):
            return kegg.name2cid(name)
    for name in name_list:
        if (kegg.name2cid(name, 0.8) != None):
            return kegg.name2cid(name, 0.8)
    return None

###########################################################################################
kegg = kegg.Kegg.getInstance()

sys.stderr.write("Parsing dG0_hatzi.csv file and adding the KEGG cids to the compounds that match by name")

csv_in = csv.reader(open('../data/dG0_hatzi.csv', 'r'))
csv_in.next()

csv_out = csv.writer(open('../rec/dG0.csv', 'w'))
for row in csv_in:
    (names, dG0, source) = row
    cid = find_cid_in_kegg(kegg, names.split('|'))
    if (cid != None):
        proper_name = kegg.cid2name(cid)
    else:
        proper_name = None
    csv_out.writerow((cid, proper_name, names, dG0, source))
