import csv
from pygibbs.thermodynamic_constants import default_c0
from toolbox.database import SqliteDatabase

class CompoundAbundance(object):
    
    def __init__(self):
        self.cids = set()
        self.media = set()
        self.default_c0 = default_c0
        self.cid2conc = {}
        self.cid2bounds = {}
    
    @staticmethod
    def StringToFloat(s):
        if not s:
            return None
        else:
            return float(s)
    
    @staticmethod
    def LoadConcentrationsFromBennett(filename=
            '../data/thermodynamics/abundance_bennett.csv'):
        abundance = CompoundAbundance()
        abundance.media = set(['glucose', 'glycerol', 'acetate'])
        for row in csv.DictReader(open(filename, 'r')):
            if not row['KEGG ID']:
                continue
            cid = int(row['KEGG ID'])
            abundance.cids.add(cid)
            if row['use'] != '1':
                continue
            
            for medium in abundance.media:
                mid_c = CompoundAbundance.StringToFloat(row[medium])
                min_c = CompoundAbundance.StringToFloat(row[medium + '(min)'])
                max_c = CompoundAbundance.StringToFloat(row[medium + '(max)'])
                if mid_c is not None:
                    abundance.cid2conc[cid, medium] = mid_c
                if min_c is not None and max_c is not None:
                    abundance.cid2bounds[cid, medium] = (min_c, max_c)
        return abundance
    
    @staticmethod
    def LoadConcentrationsFromSauer(filename=
            '../data/thermodynamics/abundance_sauer.csv'):
        """
            fields are:
            cid,glucose avg,glucose std,gluconate avg, gluconate std
            
            The concentration values are in umol/L/OD
            The conversion ratios are:
            * 3500 ul/L/OD (Volmer 2011)
            *  330 ul/L/OD (Bennett 2009)
        """
        ratio = 330 # Bennett
        
        abundance = CompoundAbundance()
        abundance.media = set(['glucose', 'gluconate'])
        for row in csv.DictReader(open(filename, 'r')):
            if not row['cid']:
                continue
            try:
                cid = int(row['cid'])
            except ValueError:
                continue
            abundance.cids.add(cid)
            for medium in abundance.media:
                avg_c = CompoundAbundance.StringToFloat(row[medium + " avg"])
                std_c = CompoundAbundance.StringToFloat(row[medium + " std"])
                if avg_c is not None:
                    abundance.cid2conc[cid, medium] = avg_c / ratio
                if avg_c is not None and std_c is not None:
                    abundance.cid2bounds[cid, medium] = \
                        ((avg_c - 2*std_c) / ratio, (avg_c + 2*std_c) / ratio)
        return abundance

    def GetAllBounds(self, medium):
        res = []
        for cid in self.cids:
            if (cid, medium) in self.cid2bounds:
                res.append([cid, self.cid2bounds[cid, medium]])
            elif (cid, medium) in self.cid2conc:
                mid_c = self.cid2conc[cid, medium]
                res.append([cid, (mid_c, mid_c)])
        return res

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