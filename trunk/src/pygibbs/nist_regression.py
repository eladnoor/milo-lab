from pygibbs.nist import Nist
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir, log_sum_exp
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import GroupDecomposer
from pygibbs import pseudoisomer
from thermodynamic_constants import R
import pylab
import logging
import csv
from toolbox.linear_regression import LinearRegression

class NistAnchors(object):
    
    def __init__(self, db, html_writer):
        self.db = db
        self.html_writer = html_writer
        self.cid2dG0_f = {}
        self.cid2min_nH = {}
    
    def FromCsvFile(self, filename='../data/thermodynamics/nist_anchors.csv'):
        self.db.CreateTable('nist_anchors', 'cid INT, charge INT, hydrogens INT, dG0 REAL')
        for row in csv.DictReader(open(filename, 'r')):
            cid = int(row['cid'])
            dG0_f = float(row['dG0'])
            z = int(row['charge'])
            nH = int(row['hydrogens'])
            self.db.Insert('nist_anchors', [cid, dG0_f, nH, z])
        
            self.cid2dG0_f[cid] = dG0_f
            self.cid2min_nH[cid] = nH
        
    def FromDatabase(self):
        for row in self.db.DictReader('nist_anchors'):
            self.cid2dG0_f[row['cid']] = row['dG0']
            self.cid2min_nH[row['cid']] = row['hydrogens']
            
    def Load(self):
        if not self.db.DoesTableExist('nist_anchors'):
            self.FromCsvFile()
        else:
            self.FromDatabase()
            
    def GetAllCids(self):
        return sorted(self.cid2dG0_f.keys())

class NistRegression(object):
    
    def __init__(self, db, html_writer, kegg):
        self.db = db
        self.html_writer = html_writer
        self.kegg = kegg
        self.nist = Nist(db, html_writer, self.kegg)
        self.nist.Load()
        
        self.nist_anchors = NistAnchors(self.db, self.html_writer)
        self.nist_anchors.FromCsvFile()
        
        dissociation = DissociationConstants(self.db, self.html_writer, self.kegg)
        dissociation.LoadValuesToDB()
        self.cid2pKa_list, self.cid2min_nH = dissociation.GetAllpKas()
        self.cid2min_charge = {}
        for cid, p in self.cid2min_nH.iteritems():
            self.cid2min_charge[cid] = self.kegg.cid2charge(cid, correctForPH=False) + \
                p - self.kegg.cid2num_hydrogens(cid, correctForPH=False)
        
    def ReverseTransform(self):
        """
            Performs the reverse Lagandre transform on all the data in NIST where
            it is possible, i.e. where all reactants have pKa values in the range
            (pH-2, pH+2) - the pH in which the Keq was measured.
        """
        cids_with_pKa = set(self.cid2pKa_list.keys())
        cids_with_pKa = set([212, 147, 9, 1, 2, 8, 13, 20])
        cids_in_nist = set(self.nist.GetAllCids())
        cids_in_nist_with_pKa = cids_with_pKa.intersection(cids_in_nist)
        
        anchored_cids = sorted(cids_in_nist_with_pKa.intersection(self.nist_anchors.GetAllCids()))
        anchored_dG0_f = pylab.matrix([[self.nist_anchors.cid2dG0_f[cid] for cid in anchored_cids]]).T
        
        unresolved_cids = sorted(cids_in_nist_with_pKa.difference(anchored_cids))
        all_cids = anchored_cids + unresolved_cids
        
        for cid in cids_with_pKa.intersection(anchored_cids):
            # TODO: Instead of raising an exception, fix the value using the pKa list
            if self.cid2min_nH[cid] != self.nist_anchors.cid2min_nH[cid]:
                raise Exception("The anchored form of C%05d is different than "
                    "the lowest form according to the pKa list %d != %d" % \
                    (cid, self.cid2min_nH[cid], self.nist_anchors.cid2min_nH[cid]))
        
        stoichiometric_matrix = pylab.zeros((0, len(all_cids)))
        dG0_r_vector = pylab.zeros((0, 1))
        
        for nist_row_data in self.nist.data:
            cids_in_reaction = set(nist_row_data.sparse.keys())
            cids_without_pKa = cids_in_reaction.difference(cids_with_pKa)
            if cids_without_pKa:
                logging.info('reaction contains CIDs with unknown pKa values: %s' % str(cids_without_pKa))
            else:
                dG0_r = nist_row_data.dG0_r + \
                    self.ReverseTransformReaction(nist_row_data.sparse, 
                    nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
                    nist_row_data.T)
                logging.debug('dG0_tag = %.1f -> dG0 = %.1f' % (nist_row_data.dG0_r, dG0_r))
                
                stoichiometric_row = pylab.zeros((1, len(all_cids)))
                for cid, coeff in nist_row_data.sparse.iteritems():
                    stoichiometric_row[0, all_cids.index(cid)] = coeff
                
                stoichiometric_matrix = pylab.vstack([stoichiometric_matrix, 
                                                      stoichiometric_row])
                dG0_r_vector = pylab.vstack([dG0_r_vector, dG0_r])
        
        
        anchored_S = stoichiometric_matrix[:, :len(anchored_cids)]
        unresolved_S = stoichiometric_matrix[:, len(anchored_cids):]
        
        unresolved_dG0_r = dG0_r_vector - anchored_S * anchored_dG0_f
        
        logging.info("Regression matrix is: %d x %d" % stoichiometric_matrix.shape)
        logging.info("%d anchored CIDs, %d unresolved CIDs" % (len(anchored_cids), len(unresolved_cids)))
        dG0_f, kerA = LinearRegression.LeastSquares(unresolved_S, unresolved_dG0_r)
        
        logging.debug("Anchors:")
        for i, cid in enumerate(anchored_cids):
            logging.debug("C%05d: dG0=%.2f, nH=%d" % (cid, anchored_dG0_f[i, 0], self.cid2min_nH[cid]))

        logging.debug("Estimations:")
        for i, cid in enumerate(unresolved_cids):
            logging.debug("C%05d: dG0=%.2f, nH=%d" % (cid, dG0_f[i, 0], self.cid2min_nH[cid]))
    
    def ReverseTransformReaction(self, sparse, pH, I, pMg, T):
        return sum([coeff * self.ReverseTransformCompound(cid, pH, I, pMg, T) \
                    for cid, coeff in sparse.iteritems()])

    def ReverseTransformCompound(self, cid, pH, I, pMg, T):
        p = self.cid2min_nH[cid]
        z = self.cid2min_charge[cid]

        # it is very important that the order of the pKa will be according
        # to increasing nH.
        pKa_list = sorted(self.cid2pKa_list[cid], reverse=True) 

        exponent_list = []
        for n in xrange(len(pKa_list)+1):
            exponent = pylab.log(10) * sum([pKa_list[i] - pH for i in xrange(n)])
            exponent += 2.91482 * ((z+n)**2 - n) * pylab.sqrt(I) / ((R * T) * (1 + 1.6 * pylab.sqrt(I)))
            exponent_list.append(exponent)
        lse = log_sum_exp(exponent_list)

        logging.debug("C%05d, p=%d, z=%d, pH=%.2f, I=%.1f, pMg=%.2f, T=%.1f, pKa=[%s], exp=[%s], lse=%.2f, correction=%.2f" % \
            (cid, p, z, pH, I, pMg, T, ','.join(['%.2f' % pKa for pKa in pKa_list]),
             ','.join(['%.2f' % pKa for pKa in exponent_list]),
             lse, R * T * (lse - p * pylab.log(10) * pH)))
        
        return R * T * (lse - p * pylab.log(10) * pH)
    
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
    logging.getLogger('').setLevel(logging.DEBUG)
    _mkdir('../res/nist/png')
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg()

    nist_regression = NistRegression(db, html_writer, kegg)
    #nist_regression.Nist_pKas()
    #nist_regression.Calculate_pKa_and_pKMg()
    
    nist_regression.ReverseTransform()
    
    html_writer.close()