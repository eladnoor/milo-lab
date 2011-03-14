import csv, logging
from kegg import Kegg
from toolbox.database import SqliteDatabase

def MatchCRCDataToKEGG():
    """
        Reads the raw data collected from the CRC handbook and tries to map every
        compound name there to a KEGG compound ID.
        Then it writes the results to the 'Public' Database
    """
    
    public_db = SqliteDatabase('../data/public_data.sqlite')
    kegg = Kegg.getInstance()

    cas2cid = {}
    for cid, comp in kegg.cid2compound_map.iteritems():
        if comp.cas:
            cas2cid[comp.cas] = cid
    
    public_db.CreateTable('pKa_from_CRC', ['cas TEXT', 'cid INT', 'name TEXT', 'formula TEXT', 'T REAL', 'pKa REAL'])
    for row_dict in csv.DictReader(open('../data/thermodynamics/pKa_from_CRC.csv', 'r')):
        cas = row_dict['CAS']
        name = row_dict['name']
        formula = row_dict['formula']
        if row_dict['T']:
            T = 273.15 + float(row_dict['T'])
        else:
            T = None
        
        if row_dict['pKa'][0] == '~':
            pKa = float(row_dict['pKa'][1:])
        else:
            pKa = float(row_dict['pKa'])

        if cas not in cas2cid:
            logging.warning('Cannot find this CAS number (%s, %s) in KEGG' % (cas, name))
            cid = None
        else:
            cid = cas2cid[cas]
            name = kegg.cid2name(cid)
        public_db.Insert('pKa_from_CRC', [unicode(cas), cid, unicode(name), unicode(formula), T, pKa])
    public_db.Commit()

if __name__ == '__main__':
    MatchCRCDataToKEGG()