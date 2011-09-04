import csv
from pygibbs.thermodynamic_constants import default_c0
from toolbox.database import SqliteDatabase

class CompoundAbundance(object):
    
    def __init__(self):
        self.media = set()
        self.default_c0 = default_c0
        self.cid2conc = {}
    
    @staticmethod
    def LoadConcentrationsFromBennett(filename=
            '../data/thermodynamics/compound_abundance.csv'):
        abundance = CompoundAbundance()
        abundance.media = set(['glucose', 'glycerol', 'acetate'])
        for row in csv.DictReader(open(filename, 'r')):
            if not row['KEGG ID']:
                continue
            cid = int(row['KEGG ID'])
            if row['use'] != '1':
                continue
            
            for medium in abundance.media:
                if row[medium] == '':
                    continue
                abundance.cid2conc[cid, medium] = float(row[medium])
        return abundance

    def ToDatabase(self, db, table_name='compound_abundance'):
        db.CreateTable(table_name, 'cid INT, medium TEXT, concentration REAL')
        for (cid, medium), conc in sorted(self.cid2conc.iteritems()):
            db.Insert(table_name, [cid, medium, conc])
        db.Commit()
    
    @staticmethod
    def FromDatabase(db, table_name='compound_abundance'):
        abundance = CompoundAbundance()
        abundance.media_list = set()
        for row in db.DictReader(table_name):
            abundance.media_list.add(row['medium'])
            abundance.cid2conc[(row['cid]'], row['medium'])] = \
                row['concentration'] # in [M]

    def GetConcentration(self, cid, c0=None, medium=None):
        c0 = c0 or self.default_c0
        if cid == 1: # the concentration of water must always be 1
            return 1
        if not medium:
            return c0 # Standard conditions = 1 [M]
        return self.cid2conc.get((cid, medium), c0)

    def GetConcentrationList(self, cid, c0=None):
        """
            return a list of pairs of media names and concentrations of the provided CID
        """
        c0 = c0 or self.default_c0
        c_list = []
        for medium in self.media_list:
            if ((cid, medium) in self.cid2conc):
                c_list.append((medium, self.cid2conc[(cid, medium)]))
        return c_list        
            
if __name__ == "__main__":
    abundance = CompoundAbundance.LoadConcentrationsFromBennett()
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    abundance.ToDatabase(db)