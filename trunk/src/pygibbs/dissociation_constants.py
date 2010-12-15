import csv, re, logging
from kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import ReadCsvWithTitles, _mkdir
from pygibbs.group_decomposition import GroupDecomposer
from toolbox.html_writer import HtmlWriter
import pybel
import sys

class DissociationConstants(object):
    def __init__(self, db, html_writer, kegg=None, group_decomposer=None):
        self.db = db
        self.html_writer = html_writer
        self.kegg = kegg or Kegg()
        self.group_decomposer = group_decomposer or GroupDecomposer.FromDatabase(db)
        self.cid2pKas = {}
    
    def ReadOriginalCSV(self, filename):
        """
            Reads the raw data collected from the CRC handbook and tries to map every
            compound name there to a KEGG compound ID.
        """
        last_formula = None
        last_name = None
        data = []
        for row_dict in ReadCsvWithTitles(filename):
            if row_dict['name'] and row_dict['formula']:
                cid, kegg_name, distance = self.FindCID(row_dict['formula'], row_dict['name'])
                last_formula = row_dict['formula']
                last_name = row_dict['name']
                row_dict['step'] = 1
            else:
                row_dict['name'] = last_name
                row_dict['formula'] = last_formula

            data.append((cid, kegg_name, distance, row_dict['name'], row_dict['formula'], row_dict['step'], row_dict['T'], row_dict['pKa']))
        
        return data
    
    def Formula2AtomBag(self, formula):
        if (formula == None or formula.find("(") != -1 or formula.find(")") != -1 or formula.find("R") != -1):
            raise Exception("non-specific compound formula: " + formula)
        
        atom_bag = {}
        for (atom, count) in re.findall("([A-Z][a-z]*)([0-9]*)", formula):
            if (count == ''):
                count = 1
            else:
                count = int(count)
            atom_bag[atom] = count

        return atom_bag
    
    def FindCID(self, formula, name):
        name = name.strip()
        cid, kegg_name, distance = self.kegg.name2cid(name, cutoff=3)
        if not cid:
            return None, None, None

        atom_bag = self.Formula2AtomBag(formula)
        kegg_atom_bag = self.kegg.cid2atom_bag(cid)
        if (atom_bag != kegg_atom_bag):
            return None, None, None
        else:
            return cid, kegg_name, distance
    
    def WriteCSV(self, fname, data):
        csv_writer = csv.writer(open(fname, 'w'))
        csv_writer.writerow(('CID', 'Kegg name', 'distance', 'original name', 'formula', 'step', 'T', 'pKa'))
        for cid, kegg_name, distance, name, formula, step, T, pKa in data:
            csv_writer.writerow((cid, kegg_name, distance, name, formula, step, T, pKa))

    def MatchCIDs(self):
        data = self.ReadOriginalCSV('../data/thermodynamics/pKa.csv')
        self.WriteCSV('../res/pKa.csv', data)
        # Now, the user must go over the output file, fix the problems in it and save it to:
        # ../data/thermodynamics/pKa_with_cids.csv
    
    def LoadValuesToDB(self, csv_filename='../data/thermodynamics/pKa_with_cids.csv'):
        """
            Load the data regarding pKa values according to KEGG compound IDs.
            First attempts to retrieve the data from the DB. If the table doesn't
            exist, it reads it from the CSV file (while caching it in the DB).
        """
        
        self.db.CreateTable('pKa', 'cid INT, step INT, T REAL, pKa REAL, smiles_below TEXT, smiles_above TEXT')
        for row in ReadCsvWithTitles(csv_filename):
            if row['cid']:
                cid = int(row['cid'])
                if row['T']:
                    T = float(row['T']) + 273.15
                else:
                    T = 298.15
                
                if row['pKa']:
                    self.db.Insert('pKa', [cid, int(row['step']), T, 
                                           float(row['pKa']), row['smiles_below'], row['smiles_above']])
                else:
                    self.db.Insert('pKa', [cid, None, None, None, None, None])
                    
        
        self.db.Commit()
            
    def AnalyseValues(self):
        for cid, step, unused_T, pKa, smiles_below, smiles_above in self.db.Execute("SELECT * FROM pKa WHERE pKa IS NOT NULL"):
            logging.info("analyzing C%05d" % cid)
            self.DrawProtonation(cid, step, pKa, smiles_below, smiles_above)
            self.cid2pKas.setdefault(cid, []).append(pKa)

        #for cid, pKas in sorted(self.cid2pKas.iteritems()):
        #    mol = self.kegg.cid2mol(cid)
        #    decomposition = self.group_decomposer.Decompose(mol, ignore_protonations=True)
        #    active_groups = self.GetActiveGroups(decomposition)
        #    print cid, self.kegg.cid2name(cid), pKas, active_groups

    def GetAllpKas(self):
        cid2pKa_list = {}
        for cid, step, unused_T, pKa, smiles_below, smiles_above in self.db.Execute("SELECT * FROM pKa"):
            cid2pKa_list.setdefault(cid, []).append(pKa)
        return cid2pKa_list   

    def GetActiveGroups(self, decomposition):
        group_name_to_index = {}

        # 'group_name_to_count' is a map from each group name to its number of appearances in 'mol'
        group_name_to_count = {}
        for i, gdata in enumerate(decomposition.groups):
            group_name, unused_protons, unused_protons, unused_mgs, node_sets = gdata
            group_name_to_index[group_name] = group_name_to_index.get(group_name, []) + [i]
            group_name_to_count[group_name] = group_name_to_count.get(group_name, 0) + len(node_sets)

        active_groups = []
        for (name, indices) in group_name_to_index.iteritems():
            if len(indices) > 1 and group_name_to_count[name] > 0:
                active_groups.append((name, group_name_to_count[name]))
        
        return active_groups
    
    def DrawProtonation(self, cid, step, pKa, smiles_below, smiles_above):
        self.html_writer.write('<h3>C%05d - %s</h3><br>\n' % (cid, self.kegg.cid2name(cid)))
        self.html_writer.write('<p>')
        if smiles_below and smiles_above:
            self.smiles2HTML(smiles_below, "C%05d_%d_below" % (cid, step))
            self.html_writer.write(" pKa = %.2f " % pKa)
            self.smiles2HTML(smiles_above, "C%05d_%d_above" % (cid, step))
        else:
            self.html_writer.write('pKa = %.2f, \n' % pKa)
            self.html_writer.write('No SMILES provided...')
        self.html_writer.write('</p>')
    
    def smiles2HTML(self, smiles, id):
        try:
            mol = pybel.readstring('smiles', str(smiles))
        except IOError:
            self.html_writer.write('Error reading smiles: ' + smiles)
            return
        mol.removeh()
        #self.html_writer.write(smiles)
        self.html_writer.embed_molecule_as_png(mol, '../res/dissociation_constants/%s.png' % id, height=50, width=50)

if (__name__ == '__main__'):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    db = SqliteDatabase("../res/gibbs.sqlite")
    html_writer = HtmlWriter("../res/dissociation_constants.html")
    _mkdir('../res/dissociation_constants')
    
    dissociation = DissociationConstants(db, html_writer)
    if False:
        dissociation.MatchCIDs()
    else:
        dissociation.LoadValuesToDB()
        dissociation.AnalyseValues()