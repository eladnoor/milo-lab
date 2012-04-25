from pygibbs.kegg import Kegg
import csv

kegg = Kegg.getInstance()
for row in csv.DictReader(open('../data/thermodynamics/abundance_bennett.csv', 'r')):
    if row['KEGG ID'] and row['reference'] == 'Rabinowitz 2009':
        cid = int(row['KEGG ID'])
        comp = kegg.cid2compound(cid)
        print ','.join([row['compound'], row['KEGG ID'], "%.1f" % comp.mass])
