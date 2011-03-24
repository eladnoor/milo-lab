import csv, logging, pybel, openbabel
from kegg import Kegg
from kegg_errors import KeggParseException
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from pygibbs.group_decomposition import GroupDecomposer
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_constants import default_T

class DissociationConstants(object):
    
    def __init__(self, db, html_writer, group_decomposer=None):
        self.db = db
        self.html_writer = html_writer
        self.kegg = Kegg.getInstance()
        self.group_decomposer = group_decomposer or GroupDecomposer.FromDatabase(db)
    
    def LoadValuesToDB(self):
        """
            Load the data regarding pKa values according to KEGG compound IDs.
        """
        
        csv_filename = '../data/thermodynamics/dissociation_constants.csv'
        self.db.CreateTable('pKa', 'cid INT, name TEXT, T REAL, nH_below INT, nH_above INT, smiles_below TEXT, smiles_above TEXT, pKa REAL')
        for line_num, row in enumerate(csv.DictReader(open(csv_filename, 'r'))):
            if row['cid']:
                cid = int(row['cid'])
            else:
                cid = None
            
            if row['T']:
                T = float(row['T'])
            else:
                T = default_T
            
            if row['pK']:
                pKa = float(row['pK'])
            else:
                pKa = None

            if row['name']:
                name = row['name']
            elif cid:
                name = self.kegg.cid2name(cid)
            else:
                name = None

            # TODO: read also the Mg dissociation constants
            if row['type'] and row['type'] != 'acid-base':
                continue
            
            try:
                self.db.Insert('pKa', [cid, name, T, int(row['nH_below']), 
                               int(row['nH_above']), 
                               row['smiles_below'], row['smiles_above'], pKa])
            except ValueError:
                raise Exception('Error while processing line #%d in %s' %\
                                (line_num, csv_filename))
                
        self.db.Commit()
            
    def AnalyseValues(self):
        for cid, nH_below, nH_above, smiles_below, smiles_above, pKa in self.db.Execute(
                "SELECT cid, nH_below, nH_above, smiles_below, smiles_above, pKa FROM pKa"
                " WHERE cid IS NOT NULL"
                " ORDER BY cid"):
            logging.info("analyzing C%05d" % cid)
            self.DrawProtonation(cid, nH_below, nH_above, smiles_below, smiles_above, pKa)

    def GetAllpKas(self):
        cid2pKa_list = {}
        cid2minimal_nH = {}
        for cid, pKa, nH_below, nH_above in self.db.Execute(
                        "SELECT cid, pKa, nH_below, nH_above FROM pKa"
                        " WHERE cid IS NOT NULL and nH_below IS NOT NULL"
                        " AND nH_above IS NOT NULL"
                        " ORDER BY cid, pKa"):
            # since we order the pKas in ascending order, the last pKa corresponds
            # to the lowest nH 
            cid2minimal_nH[cid] = nH_above
            cid2pKa_list.setdefault(cid, [])

            if pKa:
                if nH_below != nH_above+1:
                    raise Exception('The pKa=%.1f for C%05d has nH=%d below and nH=%d above it' % \
                                    (pKa, cid, nH_below, nH_above))
                cid2pKa_list[cid].append(pKa)
            elif nH_below != nH_above:
                raise Exception("When a compound has no pKa, nH_below should be"
                                " equal to nH_above, C%05d" % cid)
                
        # it is very important that the order of the pKa will be according
        # to increasing nH, i.e. decreasing pKa values
        for cid in cid2pKa_list.keys():
            cid2pKa_list[cid].sort(reverse=True)
            
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
        self.html_writer.write('<h3>C%05d - %s</h3></br>\n' % (cid, self.kegg.cid2name(cid)))
        try:
            self.html_writer.embed_molecule_as_png(self.kegg.cid2mol(cid), 'dissociation_constants/C%05d.png' % cid)
        except (KeggParseException, ValueError):
            self.html_writer.write('<b>cannot draw molecule from KEGG</b>')
        self.html_writer.write('</br>')

        if not pKa:
            self.html_writer.write('No known pKas at the physiological range')
        elif smiles_below and smiles_above:
            self.smiles2HTML(smiles_below, 'C%05d_b_H%d' % (cid, nH_below))
            self.html_writer.write(' pKa = %.2f ' % pKa)
            self.smiles2HTML(smiles_above, 'C%05d_a_H%d' % (cid, nH_above))
        elif nH_below and nH_above:
            self.html_writer.write('<b>nH = %d</b> pKa = %.2f <b>nH = %d</b>\n' % \
                                   (nH_below, pKa, nH_above))
        else:
            self.html_writer.write('<b>nH = ?</b> pKa = %.2f <b>nH = ?</b>\n' % \
                                   pKa)
        
        self.html_writer.write('</br>')
    
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
    db = SqliteDatabase("../res/gibbs.sqlite")
    html_writer = HtmlWriter("../res/dissociation_constants.html")
    _mkdir('../res/dissociation_constants')
    dissociation = DissociationConstants(db, html_writer)
    dissociation.LoadValuesToDB()
    dissociation.AnalyseValues()