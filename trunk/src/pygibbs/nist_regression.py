import logging, sys, csv
from pygibbs.nist import Nist
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import GroupDecomposer
from pygibbs import pseudoisomer


def Calculate_pKa(db, html_writer, kegg, cid2pKa_list, filename="../data/thermodynamics/dG0.csv"):
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
            pass
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

        if row['smiles']:
            smiles_dict[cid, nH, z, nMg] = row['smiles']
        else: 
            smiles_dict[cid, nH, z, nMg] = ''

    csv_writer = csv.writer(open('../res/pKa_from_dG0.csv', 'w'))
    
    html_writer.write('<table border="1">\n<tr><td>' + 
                      '</td><td>'.join(['CID', 'nH', 'charge', 'nMg', 'dG0_f', 'pKa', 'smiles_before' ,'smiles_after']) + 
                      '</td></tr>\n')
    for cid in sorted(cid2pmap.keys()):
        step = 1
        for nH, z, nMg, dG0 in sorted(cid2pmap[cid].ToMatrix(), key=lambda x:(x[2], x[0]), reverse=True):
            pKa = cid2pmap[cid].GetpKa(nH, z, nMg)
            html_writer.write('<tr><td>C%05d</td><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f</td>' % \
                              (cid, kegg.cid2name(cid), kegg.cid2formula(cid), nH, z, nMg, dG0))
            if pKa:
                html_writer.write('<td>%.2f</td><td>%s</td><td>%s</td>' \
                                  % (pKa, smiles_dict[cid, nH+1, z+1, nMg], smiles_dict[cid, nH, z, nMg]))
                if not nMg and cid not in cid2pKa_list:
                    csv_writer.writerow([cid, kegg.cid2name(cid), kegg.cid2formula(cid), step, None, "%.2f" % pKa, smiles_dict[cid, nH+1, z+1, nMg], smiles_dict[cid, nH, z, nMg]])
                    step += 1
            else:
                html_writer.write('<td>-</td><td>-</td><td>-</td>')
            html_writer.write('</tr>\n')
    html_writer.write('</table>\n')

def Nist_pKas(db, html_writer, kegg, cid2pKa_list):
    group_decomposer = GroupDecomposer.FromDatabase(db)
    cids_in_nist = set(nist.cid2count.keys())
    
    html_writer.write('CIDs with pKa: %d<br>\n' % len(cid2pKa_list))
    html_writer.write('CIDs in NIST: %d<br>\n' % len(cids_in_nist))
    html_writer.write('CIDs in NIST with pKas: %d<br>\n' % len(cids_in_nist.intersection(cid2pKa_list.keys())))
    
    html_writer.write('All CIDs in NIST: <br>\n')
    html_writer.write('<table border="1">\n')
    html_writer.write('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' % ("CID", "NAME", "COUNT", "REMARK"))
    for cid, count in sorted(nist.cid2count.iteritems()):
        if cid not in cid2pKa_list:
            html_writer.write('<tr><td>C%05d</td><td>%s</td><td>%d</td><td>' % (cid, kegg.cid2name(cid), count))
            try:
                mol = kegg.cid2mol(cid)
                decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
    
                if len(decomposition.PseudoisomerVectors()) > 1:
                    html_writer.write('should have pKas')
                else:
                    html_writer.write('doesn\'t have pKas')
                html_writer.embed_molecule_as_png(kegg.cid2mol(cid), 'png/C%05d.png' % cid)
            
            except Exception:
                html_writer.write('cannot decompose')
            html_writer.write('</td></tr>\n')
    
    html_writer.write('</table>\n')

def get_cid2pKa_list(db, html_writer, kegg):
    dissociation = DissociationConstants(db, html_writer, kegg)
    dissociation.LoadValuesToDB('../data/thermodynamics/pKa_with_cids.csv')
    return dissociation.GetAllpKas()

if (__name__ == "__main__"):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    _mkdir('../res/nist/png')
    
    kegg = Kegg()
    nist = Nist(kegg)
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    
    cid2pKa_list = get_cid2pKa_list(db, html_writer, kegg)
    Nist_pKas(db, html_writer, kegg, cid2pKa_list)
    Calculate_pKa(db, html_writer, kegg, cid2pKa_list)
    
    html_writer.close()