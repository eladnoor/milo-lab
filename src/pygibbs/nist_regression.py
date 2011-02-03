#!/usr/bin/python

import pylab
import logging
import csv
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist import Nist
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecomposer
from pygibbs.thermodynamics import Thermodynamics,\
    MissingCompoundFormationEnergy, CsvFileThermodynamics
from pygibbs.pseudoisomers_data import DissociationTable
from pygibbs.pseudoisomer import PseudoisomerMap

class NistAnchors(object):
    
    def __init__(self, db, html_writer):
        self.db = db
        self.html_writer = html_writer
        self.cid2dG0_f = {}
        self.cid2min_nH = {}
        
    def FromCsvFile(self, filename='../data/thermodynamics/nist_anchors.csv'):
        self.db.CreateTable('nist_anchors', 'cid INT, z INT, nH INT, nMg INT, dG0 REAL')
        for row in csv.DictReader(open(filename, 'r')):
            cid = int(row['cid'])
            dG0_f = float(row['dG0'])
            z = int(row['z'])
            nH = int(row['nH'])
            nMg = int(row['nMg'])
            self.db.Insert('nist_anchors', [cid, dG0_f, nH, z, nMg])
        
            self.cid2dG0_f[cid] = dG0_f
            self.cid2min_nH[cid] = nH
        
    def FromDatabase(self):
        for row in self.db.DictReader('nist_anchors'):
            self.cid2dG0_f[row['cid']] = row['dG0']
            self.cid2min_nH[row['cid']] = row['nH']
            
    def Load(self):
        if not self.db.DoesTableExist('nist_anchors'):
            self.FromCsvFile()
        else:
            self.FromDatabase()
            
    def GetAllCids(self):
        return sorted(self.cid2dG0_f.keys())

class NistRegression(Thermodynamics):
    
    def __init__(self, db, html_writer, kegg=None, nist=None):
        Thermodynamics.__init__(self)
        self.db = db
        self.html_writer = html_writer
        self.kegg = kegg or Kegg()
        if nist:
            self.nist = nist
        else:
            self.nist = Nist(db, html_writer, self.kegg)
            self.nist.Load()
        
        self.nist_anchors = NistAnchors(self.db, self.html_writer)
        self.nist_anchors.FromCsvFile()
        
        self.cid2diss_table = DissociationTable.ReadDissociationCsv(kegg=kegg)
        self.cid2pmap_dict = {}
        
    def cid2PseudoisomerMap(self, cid):
        if (cid in self.cid2pmap_dict):
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())
        
    def ReverseTransform(self, prior_thermodynamics=None, T_range=None, assume_no_pKa_by_default=False):
        """
            Performs the reverse Lagandre transform on all the data in NIST where
            it is possible, i.e. where all reactants have pKa values in the range
            (pH-2, pH+2) - the pH in which the Keq was measured.
        """
        cids_with_pKa = set(self.cid2diss_table.keys())
        cids_in_nist = set(self.nist.GetAllCids())
        
        if assume_no_pKa_by_default:
            for cid in cids_in_nist.difference(cids_with_pKa):
                diss = DissociationTable(cid)
                diss.min_nH = self.kegg.cid2num_hydrogens(cid)
                diss.min_charge = self.kegg.cid2charge(cid)
                if diss.min_nH == None or diss.min_charge == None:
                    logging.warning('cannot add C%05d since nH or charge '
                                    'cannot be determined' % cid) 
                else:
                    self.cid2diss_table[cid] = diss
                    cids_with_pKa.add(cid)
        
        cids_in_nist_with_pKa = cids_with_pKa.intersection(cids_in_nist)
        logging.info('There are pKa values for %d out of the %d compounds in NIST' % \
                     (len(cids_in_nist_with_pKa), len(cids_in_nist)))

        #self.anchors = set(cids_in_nist_with_pKa.intersection(self.nist_anchors.GetAllCids()))
        cids_to_estimate = sorted(cids_in_nist_with_pKa.difference(self.anchors))
        logging.info("%d compounds are anchored, and the %d others will be estimated" % \
                     (len(self.anchors), len(cids_to_estimate)))
        
        # get a vector of anchored formation energies. one needs to be careful
        # to always use the most basic pseudoisomer (the one with the lowest nH)
        # because these are the forms used in the regression matrix
        for cid in self.anchors:
            dG0_f = self.nist_anchors.cid2dG0_f[cid]
            nH = self.nist_anchors.cid2min_nH[cid]
            dG0_f_base = self.ConvertPseudoisomer(cid, dG0_f, nH)

            # by assigning a dG0_f to this compound, it will be subtracted out
            # from the dG0_r when applying the reverse transform.
            self.cid2diss_table[cid].min_dG0 = dG0_f_base
            self.cid2pmap_dict[cid] = self.cid2diss_table[cid].GetPseudoisomerMap()

        stoichiometric_matrix = pylab.zeros((0, len(cids_to_estimate)))
        
        dG0_r_tag = pylab.zeros((0, 1))
        ddG0_r = pylab.zeros((0, 1)) # the difference between dG0_r and dG'0_r
        
        nist_rows_used = []
        for nist_row_data in self.nist.data:
            if T_range and not (T_range[0] < nist_row_data.T < T_range[1]):
                logging.debug('Temperature %.2f not within allowed range', nist_row_data.T)
                continue # the temperature is outside the allowed range
            
            cids_in_reaction = set(nist_row_data.sparse.keys())
            cids_without_pKa = cids_in_reaction.difference(cids_with_pKa)
            if cids_without_pKa:
                logging.debug('reaction contains CIDs with unknown pKa values: %s' % \
                              ', '.join(['C%05d' % cid for cid in cids_without_pKa]))
            else:
                nist_rows_used.append(nist_row_data)
                dG0_r_tag = pylab.vstack([dG0_r_tag, nist_row_data.dG0_r])
                
                ddG = self.ReverseTransformReaction(nist_row_data.sparse, 
                    nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
                    nist_row_data.T)
                ddG0_r = pylab.vstack([ddG0_r, ddG])
                
                stoichiometric_row = pylab.zeros((1, len(cids_to_estimate)))
                for cid, coeff in nist_row_data.sparse.iteritems():
                    if cid not in self.anchors:
                        stoichiometric_row[0, cids_to_estimate.index(cid)] = coeff
                
                stoichiometric_matrix = pylab.vstack([stoichiometric_matrix, 
                                                      stoichiometric_row])
        
        logging.info("%d out of %d NIST measurements can be used" % \
                     (stoichiometric_matrix.shape[0], len(self.nist.data)))

        dG0_r = dG0_r_tag + ddG0_r
        # squeeze the regression matrix by leaving only unique rows
        unique_rows_S = pylab.unique([tuple(stoichiometric_matrix[i,:].flat) for i 
                                      in xrange(stoichiometric_matrix.shape[0])])

        logging.info("There are %d unique reactions" % \
                     unique_rows_S.shape[0])
        
        # for every unique row, calculate the average dG0_r of all the rows that
        # are the same reaction
        unique_rows_dG0_r = pylab.zeros((0, 1))
        dG0_r = pylab.hstack([pylab.zeros(dG0_r.shape), dG0_r])
        for i in xrange(unique_rows_S.shape[0]):
            # find the list of indices which are equal to row i in unique_rows_S
            diff = abs(stoichiometric_matrix - unique_rows_S[i,:])
            row_indices = pylab.find(pylab.sum(diff, 1)==0)
            
            # take the average of the dG0_r of these rows
            average_dG0_r = pylab.mean([dG0_r[j, 1] for j in row_indices])
            unique_rows_dG0_r = pylab.vstack([unique_rows_dG0_r, average_dG0_r])
            dG0_r[row_indices, 0] = average_dG0_r
            
        fig = pylab.figure()
        pylab.plot(dG0_r[:,0], dG0_r[:,1]-dG0_r[:,0], '.')
        pylab.xlabel("$<\Delta_r G^\circ>$")
        pylab.ylabel("$\Delta_r G^\circ - <\Delta_r G^\circ>$")
        pylab.title('$\sigma = %.1f$ kJ/mol' % (pylab.std(dG0_r[:,1]-dG0_r[:,0])))
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480)

        pylab.np.savetxt('../res/nist/regress_CID.txt', 
            pylab.array(cids_to_estimate), fmt='%d', delimiter=',')
        pylab.np.savetxt('../res/nist/regress_S.txt', 
            unique_rows_S, fmt='%g', delimiter=',')
        pylab.np.savetxt('../res/nist/regress_dG0.txt',
            unique_rows_dG0_r, fmt='%.2f', delimiter=',')
        
        logging.info("Regression matrix is %d x %d, and it's rank is %d" % \
                     (unique_rows_S.shape[0], unique_rows_S.shape[1],
                      LinearRegression.Rank(unique_rows_S)))
        estimated_dG0_f, kerA = LinearRegression.LeastSquares(unique_rows_S, 
            unique_rows_dG0_r, reduced_row_echlon=False)
        corr = pylab.corrcoef(pylab.dot(unique_rows_S, estimated_dG0_f), 
                              unique_rows_dG0_r, rowvar=0)[0, 1]
        logging.info("Regression Complete, r^2 = %.6f" % corr**2)
        logging.info("The dimension of the Kernel is %d" % (kerA.shape[0]))

        if prior_thermodynamics:
            # find the vector in the solution subspace which is closest to the 
            # prior formation energies
            delta_dG0_f = pylab.zeros((0, 1))
            indices_in_prior = []
            for i, cid in enumerate(cids_to_estimate):
                try:
                    pmap = prior_thermodynamics.cid2PseudoisomerMap(cid)
                    for p_nH, unused_z, p_nMg, dG0 in sorted(pmap.ToMatrix()):
                        if p_nMg == 0:
                            dG0_base = self.ConvertPseudoisomer(cid, dG0, p_nH)
                            difference = dG0_base - estimated_dG0_f[i, 0]
                            delta_dG0_f = pylab.vstack([delta_dG0_f, difference])
                            indices_in_prior.append(i)
                            break
                except MissingCompoundFormationEnergy:
                    continue
            
            v, _ = LinearRegression.LeastSquares(kerA.T[indices_in_prior,:], 
                        delta_dG0_f, reduced_row_echlon=False)
            estimated_dG0_f += pylab.dot(kerA.T, v)

        # copy the solution into the diss_tables of all the compounds,
        # and then generate their PseudoisomerMaps.
        for i, cid in enumerate(cids_to_estimate):
            self.cid2diss_table[cid].min_dG0 = estimated_dG0_f[i, 0]
            self.cid2pmap_dict[cid] = self.cid2diss_table[cid].GetPseudoisomerMap()
            
    def ConvertPseudoisomer(self, cid, dG0, nH_from, nH_to=None):
        try:
            return self.cid2diss_table[cid].ConvertPseudoisomer(dG0, nH_from, nH_to)
        except KeyError as e:
            raise KeyError("In C%05d, %s" % (cid, str(e)))
    
    def ReverseTransformReaction(self, sparse, pH, I, pMg, T):
        return sum([coeff * self.ReverseTransformCompound(cid, pH, I, pMg, T) \
                    for cid, coeff in sparse.iteritems()])

    def ReverseTransformCompound(self, cid, pH, I, pMg, T):
        return -self.cid2diss_table[cid].Transform(pH, I, pMg, T)

    def Nist_pKas(self):
        group_decomposer = GroupDecomposer.FromDatabase(self.db)
        cids_in_nist = set(self.nist.cid2count.keys())
        cids_with_pKa = set(self.cid2diss_table.keys())
        
        self.html_writer.write('CIDs with pKa: %d<br>\n' % len(cids_with_pKa))
        self.html_writer.write('CIDs in NIST: %d<br>\n' % len(cids_in_nist))
        self.html_writer.write('CIDs in NIST with pKas: %d<br>\n' % \
                          len(cids_in_nist.intersection(cids_with_pKa)))
        
        self.html_writer.write('All CIDs in NIST: <br>\n')
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' % ("cid", "name", "count", "remark"))
        for cid, count in sorted(self.nist.cid2count.iteritems()):
            if cid not in cids_with_pKa:
                self.html_writer.write('<tr><td><a href="%s">C%05d<a></td><td>%s</td><td>%d</td><td>' % \
                    (self.kegg.cid2link(cid), cid, self.kegg.cid2name(cid), count))
                try:
                    mol = self.kegg.cid2mol(cid)
                    decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
        
                    if len(decomposition.PseudoisomerVectors()) > 1:
                        self.html_writer.write('should have pKas')
                    else:
                        self.html_writer.write('doesn\'t have pKas')
                    self.html_writer.embed_molecule_as_png(
                        self.kegg.cid2mol(cid), 'png/C%05d.png' % cid)
                
                except Exception:
                    self.html_writer.write('cannot decompose')
                self.html_writer.write('</td></tr>\n')
        
        self.html_writer.write('</table>\n')

    def Calculate_pKa_and_pKMg(self, filename="../data/thermodynamics/dG0.csv"):
        cid2pmap = {}
        smiles_dict = {}
        
        for row in csv.DictReader(open(filename, 'r')):
            #smiles, cid, compound_name, dG0, unused_dH0, charge, hydrogens, Mg, use_for, ref, unused_assumption 
            name = "%s (z=%s, nH=%s, nMg=%s)" % (row['compound name'], row['z'], row['nH'], row['nMg'])
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
                    nH = int(row['nH'])
                    z = int(row['z'])
                    nMg = int(row['nMg'])
                except ValueError:
                    raise Exception("can't read the data about %s" % (row['compound name']))
                cid2pmap.setdefault(cid, PseudoisomerMap())
                cid2pmap[cid].Add(nH, z, nMg, dG0)
    
            if row['smiles']:
                smiles_dict[cid, nH, z, nMg] = row['smiles']
            else: 
                smiles_dict[cid, nH, z, nMg] = ''
    
        #csv_writer = csv.writer(open('../res/pKa_from_dG0.csv', 'w'))
        
        self.self.html_writer.write('<table border="1">\n<tr><td>' + 
                          '</td><td>'.join(['cid', 'name', 'formula', 'nH', 'z', 'nMg', 'dG0_f', 'pKa', 'pK_Mg']) + 
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

    def ToDatabase(self):
        Thermodynamics.ToDatabase(self, self.db, 'nist_regression')

    def FromDatabase(self):
        if self.db.DoesTableExist('nist_regression'):
            Thermodynamics.FromDatabase(self, self.db, 'nist_regression')
        else:
            logging.warning('You should run nist_regression.py before trying to'
                            ' load the data from the database')
        
    def WriteDataToHtml(self):
        Thermodynamics.WriteDataToHtml(self, self.html_writer, self.kegg)
        
    def VerifyResults(self, T_range=None):
        return self.nist.verify_results(self, T_range)


def main():
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    alberty = CsvFileThermodynamics('../data/thermodynamics/alberty_pseudoisomers.csv')
    alberty.ToDatabase(db, 'alberty')
    
    html_writer.write("<h2>NIST regression:</h2>")
    nist_regression = NistRegression(db, html_writer, kegg)
    
    if False:
        nist_regression.Nist_pKas()
        #nist_regression.Calculate_pKa_and_pKMg()
    else:
        T_range = (298, 314)
        
        nist_regression.ReverseTransform(prior_thermodynamics=alberty,
                                         T_range=T_range,
                                         assume_no_pKa_by_default=True)

        nist_regression.ToDatabase()
        
        html_writer.write('<h3>Regression results:</h3>\n')
        html_writer.insert_toggle('regression')
        html_writer.start_div('regression')
        nist_regression.WriteDataToHtml()
        html_writer.end_div()
    
        html_writer.write('<h3>Reaction energies - Estimated vs. Observed:</h3>\n')
        html_writer.insert_toggle('verify')
        html_writer.start_div('verify')
        N, rmse = nist_regression.VerifyResults(T_range)
        html_writer.end_div()
        html_writer.write('</br>\n')
        
        logging.info("N = %d, RMSE = %.1f" % (N, rmse))

        html_writer.write('<h3>Formation energies - Estimated vs. Alberty:</h3>\n')

        query = 'SELECT k.cid, k.name, a.nH, a.z, a.nMg, a.dG0_f, r.dG0_f ' + \
                'FROM kegg_compound k, alberty a, nist_regression r ' + \
                'WHERE k.cid=a.cid AND a.cid=r.cid AND a.nH=r.nH AND a.nMg=r.nMg ' + \
                'AND a.anchor=0 ORDER BY a.cid,a.nH'
        column_names = ['CID', 'name', 'nH', 'z', 'nMg', 'dG0_f(Alberty)',
                        'dG0_f(Regression)']
        
        data = pylab.zeros((0, 2))
        fig = pylab.figure()
        pylab.hold(True)
        for row in db.Execute(query):
            unused_cid, name, unused_nH, z, unused_nMg, dG0_a, dG0_r = row
            x = (dG0_a + dG0_r)/2
            y = dG0_a - dG0_r
            pylab.text(x, y, "%s [%d]" % (name, z), fontsize=5, rotation=20)
            data = pylab.vstack([data, (x,y)])

        pylab.plot(data[:,0], data[:,1], '.')
        html_writer.embed_matplotlib_figure(fig, width=640, height=480)
        db.Query2CSV('../res/nist/alberty_vs_regress.csv', query, column_names)

    html_writer.close()
    
if (__name__ == "__main__"):
    main()