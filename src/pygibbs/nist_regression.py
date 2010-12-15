import logging, sys, csv
from pygibbs.nist import Nist
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import GroupDecomposer
from pygibbs import pseudoisomer


def Calculate_pKa(db, html_writer, kegg, cids_with_pKa, filename="../data/thermodynamics/dG0.csv"):
    cid2pmap = {}
    smiles_dict = {}
    
    for row in csv.DictReader(open(filename, 'r')):
        #smiles, cid, compound_name, dG0, unused_dH0, charge, hydrogens, Mg, use_for, ref, unused_assumption 
        name = "%s (z=%s, nH=%s, nMg=%s)" % (row['compound name'], row['charge'], row['hydrogens'], row['Mg'])
        logging.info('reading data for ' + name)

        if not row['dG0']:
            continue

        if (row['use for'] == "skip"):
            continue
            
        try:
            dG0 = float(row['dG0'])
        except ValueError:
            raise Exception("Invalid dG0: " + str(dG0))

        if (row['use for'] == "test"):
            continue
        elif (row['use for'] == "train"):
            pass
        else:
            raise Exception("Unknown usage flag: " + row['use for'])

        if row['cid']:
            cid = int(row['cid'])
            try:
                nH = int(row['hydrogens'])
                z = int(row['charge'])
                nMg = int(row['Mg'])
            except ValueError:
                raise Exception("can't read the data about %s" % (row['compound name']))
            cid2pmap.setdefault(cid, pseudoisomer.PseudoisomerMap())
            cid2pmap[cid].Add(nH, z, nMg, dG0)

        if (row['smiles'] == ""):
            raise Exception("Cannot use compound '%s' for training if it lacks a SMILES string" % row['compound name'])
        smiles_dict[cid, nH, z, nMg] = row['smiles'] 

    csv_writer = csv.writer(open('../res/pKa_from_dG0.csv', 'w'))
    
    html_writer.write('<table border="1">\n<tr><td>' + 
                      '</td><td>'.join(['CID', 'nH', 'charge', 'nMg', 'dG0_f', 'pKa', 'smiles_before' ,'smiles_after']) + 
                      '</td></tr>\n')
    for cid in sorted(cid2pmap.keys()):
        step = 1
        for nH, z, nMg, dG0 in sorted(cid2pmap[cid].ToMatrix(), key=lambda x:(x[2], x[0]), reverse=True):
            pKa = cid2pmap[cid].GetpKa(nH, z, nMg)
            if pKa:
                html_writer.write('<tr><td>C%05d</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f</td><td>%.2f</td><td>%s</td><td>%s</td></tr>\n' \
                                  % (cid, nH, z, nMg, dG0, pKa, smiles_dict[cid, nH+1, z+1, nMg], smiles_dict[cid, nH, z, nMg]))
                if not nMg and cid not in cids_with_pKa:
                    csv_writer.writerow([cid, kegg.cid2name(cid), kegg.cid2formula(cid), step, None, "%.2f" % pKa, smiles_dict[cid, nH+1, z+1, nMg], smiles_dict[cid, nH, z, nMg]])
                    step += 1
            else:
                html_writer.write('<tr><td>C%05d</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f</td><td>-</td><td>-</td><td>-</td></tr>\n' \
                                  % (cid, nH, z, nMg, dG0))
    html_writer.write('</table>\n')

def Nist_pKas(db, html_writer):
    group_decomposer = GroupDecomposer.FromDatabase(db)
    dissociation = DissociationConstants(db, html_writer, kegg)
    dissociation.LoadValuesToDB('../data/thermodynamics/pKa_with_cids.csv')
    cid2pKa_list = dissociation.GetAllpKas()
    
    cids_with_pKa = set(cid2pKa_list.keys())
    cids_in_nist = set(nist.cid2count.keys())
    
    html_writer.write("CIDs with pKa: %d<br>\n" % len(cids_with_pKa))
    html_writer.write("CIDs in NIST: %d<br>\n" % len(cids_in_nist))
    html_writer.write("CIDs in NIST with pKas: %d<br>\n" % len(cids_in_nist.intersection(cids_with_pKa)))
    
    html_writer.write("All CIDs in NIST: <br>\n<ul>\n")
    for cid in sorted(cids_in_nist):
        if cid in cids_with_pKa:
            html_writer.write("  <li>C%05d - has pKas</li>" % cid)
            continue
        try:
            mol = kegg.cid2mol(cid)
            decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
        except Exception:
            html_writer.write("  <li>C%05d - cannot decompose</li>" % cid)
            continue
        
        if len(decomposition.PseudoisomerVectors()) > 1:
            html_writer.write("  <li>C%05d - should have pKas</li>" % cid)
        else:
            html_writer.write("  <li>C%05d - doesn't have pKas</li>" % cid)
    html_writer.write("</ul>\n")
    return cids_with_pKa

if (__name__ == "__main__"):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    _mkdir("../res/nist")
    
    kegg = Kegg()
    nist = Nist(kegg)
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    
    cids_with_pKa = Nist_pKas(db, html_writer)
    Calculate_pKa(db, html_writer, kegg, cids_with_pKa)
    
    html_writer.close()