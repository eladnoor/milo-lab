from pygibbs.nist import Nist
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import GroupDecomposer
from pygibbs import pseudoisomer
from thermodynamic_constants import R
import pylab
import logging
import csv

class NistRegression(object):
    
    def __init__(self, db, html_writer, kegg):
        self.db = db
        self.html_writer = html_writer
        self.kegg = kegg
        self.nist = Nist(db, html_writer, self.kegg)
        self.nist.FromDatabase()
        dissociation = DissociationConstants(self.db, self.html_writer, self.kegg)
        dissociation.LoadValuesToDB('../data/thermodynamics/pKa_with_cids.csv')
        self.cid2pKa_list, self.cid2min_nH = dissociation.GetAllpKas()
        
    def ReverseTransform(self):
        """
            Performs the reverse Lagandre transform on all the data in NIST where
            it is possible, i.e. where all reactants have pKa values in the range
            (pH-2, pH+2) - the pH in which the Keq was measured.
        """
        for nist_row_data in self.nist.data():
            dG0 = nist_row_data.dG0 + \
                self.ReverseTransformReaction(nist_row_data.sparse, 
                nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
                nist_row_data.T)
            logging.info('dG0_tag = %.1f -> dG0 = %.1f' % (nist_row_data.dG0, dG0))
    
    def ReverseTransformReaction(self, sparse, pH, I, pMg, T):
        return sum([coeff * self.ReverseTransformCompound(cid, pH, I, pMg, T) \
                    for coeff, cid in sparse.iteritems()])

    def ReverseTransformCompound(self, cid, pH, I, pMg, T):
        sum_exp = 0
        pKa_list = self.cid2pKa_list[cid]
        for n in xrange(1, len(pKa_list)):
            sum_exp += 10**sum([pH - pKa_list[i] for i in xrange(n)])

        p = self.cid2min_nH[cid]

        return R * T * (pylab.log(sum_exp) - p * pylab.log(10) * pH)
    
    def Nist_pKas(self):
        group_decomposer = GroupDecomposer.FromDatabase(self.db)
        cids_in_nist = set(self.nist.cid2count.keys())
        
        self.html_writer.write('CIDs with pKa: %d<br>\n' % len(self.cid2pKa_list))
        self.html_writer.write('CIDs in NIST: %d<br>\n' % len(cids_in_nist))
        self.html_writer.write('CIDs in NIST with pKas: %d<br>\n' % \
                          len(cids_in_nist.intersection(self.cid2pKa_list.keys())))
        
        self.html_writer.write('All CIDs in NIST: <br>\n')
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' % ("CID", "NAME", "COUNT", "REMARK"))
        for cid, count in sorted(self.nist.cid2count.iteritems()):
            if cid not in self.cid2pKa_list:
                self.self.html_writer.write('<tr><td>C%05d</td><td>%s</td><td>%d</td><td>' % (cid, self.kegg.cid2name(cid), count))
                try:
                    mol = self.kegg.cid2mol(cid)
                    decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
        
                    if len(decomposition.PseudoisomerVectors()) > 1:
                        self.self.html_writer.write('should have pKas')
                    else:
                        self.self.html_writer.write('doesn\'t have pKas')
                    self.self.html_writer.embed_molecule_as_png(
                        self.kegg.cid2mol(cid), 'png/C%05d.png' % cid)
                
                except Exception:
                    self.self.html_writer.write('cannot decompose')
                self.self.html_writer.write('</td></tr>\n')
        
        self.self.html_writer.write('</table>\n')

    def Calculate_pKa_and_pKMg(self, filename="../data/thermodynamics/dG0.csv"):
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
    
        #csv_writer = csv.writer(open('../res/pKa_from_dG0.csv', 'w'))
        
        self.self.html_writer.write('<table border="1">\n<tr><td>' + 
                          '</td><td>'.join(['CID', 'name', 'formula', 'nH', 'charge', 'nMg', 'dG0_f', 'pKa', 'pK_Mg']) + 
                          '</td></tr>\n')
        for cid in sorted(cid2pmap.keys()):
            #step = 1
            for nH, z, nMg, dG0 in sorted(cid2pmap[cid].ToMatrix(), key=lambda x:(-x[2], -x[0])):
                pKa = cid2pmap[cid].GetpKa(nH, z, nMg)
                pK_Mg = cid2pmap[cid].GetpK_Mg(nH, z, nMg)
                self.self.html_writer.write('<tr><td>')
                self.self.html_writer.write('</td><td>'.join(["C%05d" % cid, 
                    self.kegg.cid2name(cid) or "?", 
                    self.kegg.cid2formula(cid) or "?", 
                    str(nH), str(z), str(nMg), 
                    "%.1f" % dG0, str(pKa), str(pK_Mg)]))
                #if not nMg and cid not in cid2pKa_list:
                #    csv_writer.writerow([cid, kegg.cid2name(cid), kegg.cid2formula(cid), step, None, "%.2f" % pKa, smiles_dict[cid, nH+1, z+1, nMg], smiles_dict[cid, nH, z, nMg]])
                #    step += 1
                self.self.html_writer.write('</td></tr>\n')
        self.self.html_writer.write('</table>\n')

if (__name__ == "__main__"):
    _mkdir('../res/nist/png')
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg()

    nist_regression = NistRegression(db, html_writer, kegg)
    #nist_regression.Nist_pKas()
    #nist_regression.Calculate_pKa_and_pKMg()
    
    nist_regression.ReverseTransform()
    
    html_writer.close()