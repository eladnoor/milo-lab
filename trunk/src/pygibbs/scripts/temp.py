import csv
from toolbox.molecule import Molecule


filename = '../data/thermodynamics/dissociation_constants.csv'
csv_reader = csv.DictReader(open(filename, 'r'))

csv_writer = csv.writer(open('../res/temp.csv', 'w'))
cid = 90000
name2cid = {}
for i, row in enumerate(csv_reader):
    if not row['cid']:
        name = row['name']
        if name in name2cid:
            continue
        cid += 1
        name2cid[name] = cid
        nH_below = int(row['nH_below'])
        nH_above = int(row['nH_above'])
        nMg_below = int(row['nMg_below'])
        nMg_above = int(row['nMg_above'])
        smiles_below = row['smiles_below']
        smiles_above = row['smiles_above']

        mol = Molecule.FromSmiles(smiles_below)
        inchi = mol.ToInChI()
        
        csv_writer.writerow([name, cid, inchi]) 
        