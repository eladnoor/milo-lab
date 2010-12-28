import csv, re, logging
from kegg import Kegg, KeggParseException
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from pygibbs.group_decomposition import GroupDecomposer
from toolbox.html_writer import HtmlWriter
import pybel
import sys
import openbabel
from pygibbs.thermodynamic_constants import default_T

class DissociationConstants(object):
    def __init__(self, db, html_writer, kegg=None, group_decomposer=None):
        self.db = db
        self.html_writer = html_writer
        self.kegg = kegg or Kegg()
        self.group_decomposer = group_decomposer or GroupDecomposer.FromDatabase(db)
    
    def ReadOriginalCSV(self, filename):
        """
            Reads the raw data collected from the CRC handbook and tries to map every
            compound name there to a KEGG compound ID.
        """
        last_formula = None
        last_name = None
        data = []
        for row_dict in csv.DictReader(open(filename, 'r')):
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
    
    def LoadValuesToDB(self, csv_filename='../data/thermodynamics/pKa_with_nH.csv'):
        """
            Load the data regarding pKa values according to KEGG compound IDs.
        """
        
        self.db.CreateTable('pKa', 'cid INT, T REAL, nH_below INT, nH_above INT, smiles_below TEXT, smiles_above TEXT, pKa REAL')
        for row in csv.DictReader(open(csv_filename, 'r')):
            cid = int(row['CID'])
            if row['T']:
                T = float(row['T']) + 273.15
            else:
                T = default_T
            
            pKa = row['pKa'] or None
            if pKa:
                pKa = float(pKa)
            nH_below = row['nH_below'] or None
            nH_above = row['nH_above'] or None
            self.db.Insert('pKa', [cid, T, nH_below, nH_above, 
                                   row['smiles_below'], row['smiles_above'],
                                   pKa])
                
        self.db.Commit()
            
    def AnalyseValues(self):
        #csv_writer = csv.writer(open('../res/pKa_with_nH.csv', 'w'))
        #csv_writer.writerow(['CID', 'T', 'nH_below', 'nH_above', 'smiles_below', 'smiles_above', 'pKa'])
        for cid, nH_below, nH_above, smiles_below, smiles_above, pKa in self.db.Execute(
                "SELECT cid, nH_below, nH_above, smiles_below, smiles_above, pKa FROM pKa ORDER BY cid"):
            logging.info("analyzing C%05d" % cid)
            self.DrawProtonation(cid, nH_below, nH_above, smiles_below, smiles_above, pKa)
            #if smiles_below and smiles_above:
            #    nH_below = DissociationConstants.smiles2nH(smiles_below)
            #    nH_above = DissociationConstants.smiles2nH(smiles_above)
            #else:
            #    nH_below = None
            #    nH_above = self.kegg.cid2num_hydrogens(cid)
            #csv_writer.writerow([cid, T, nH_below, nH_above, smiles_below, smiles_above, pKa])

    def GetAllpKas(self):
        cid2pKa_list = {}
        cid2minimal_nH = {}
        for cid, pKa, nH_below, nH_above in self.db.Execute(
                        "SELECT cid, pKa, nH_below, nH_above FROM pKa"):
            cid2pKa_list.setdefault(cid, [])
            if pKa:
                if nH_below != nH_above+1:
                    raise Exception('The pKa=%.1f for C%05d has nH=%d below and nH=%d above it' % \
                                    (pKa, cid, nH_below, nH_above))
                cid2pKa_list[cid].append(pKa)
                if pKa == max(cid2pKa_list[cid]):
                    cid2minimal_nH[cid] = nH_above
            else:
                cid2minimal_nH[cid] = self.kegg.cid2num_hydrogens(cid)
                
        return cid2pKa_list, cid2minimal_nH

    @staticmethod
    def smiles2nH(smiles, correctForPH=False):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smiles", "mol")
        obmol = openbabel.OBMol()
        polaronly = False
        obConversion.ReadString(obmol, str(smiles))
        obmol.AddHydrogens(polaronly, correctForPH)
        return obmol.NumAtoms() - obmol.NumHvyAtoms()

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
    
    def DrawProtonation(self, cid, nH_below, nH_above, smiles_below, smiles_above, pKa):
        self.html_writer.write('<h3>C%05d - %s</h3><br>\n' % (cid, self.kegg.cid2name(cid)))
        self.html_writer.write('<p>')
        if not pKa:
            self.html_writer.write('No known pKas at the physiological range')
            self.html_writer.embed_molecule_as_png(self.kegg.cid2mol(cid), 'dissociation_constants/C%05d.png' % cid)
        elif smiles_below and smiles_above:
            self.smiles2HTML(smiles_below, "C%05d_b_H%d" % (cid, nH_below))
            self.html_writer.write(" pKa = %.2f " % pKa)
            self.smiles2HTML(smiles_above, "C%05d_a_H%d" % (cid, nH_above))
        else:
            self.html_writer.write('pKa = %.2f, \n' % pKa)
            self.html_writer.write('no SMILES provided...')
            try:
                self.html_writer.embed_molecule_as_png(self.kegg.cid2mol(cid), 'dissociation_constants/C%05d.png' % cid)
            except KeggParseException:
                pass
        self.html_writer.write('</p>')
    
    def smiles2HTML(self, smiles, id, height=50, width=50):
        try:
            mol = pybel.readstring('smiles', str(smiles))
        except IOError:
            self.html_writer.write('Error reading smiles: ' + smiles)
            return
        mol.removeh()
        #self.html_writer.write(smiles)
        self.html_writer.embed_molecule_as_png(mol, 'dissociation_constants/%s.png' % id, height=height, width=width)

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