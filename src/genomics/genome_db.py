import csv
import json

from toolbox import database
from toolbox import util

from pygibbs.kegg import Kegg


class GenomeDB(object):
    
    OXY_REQ = 'Oxygen Requirement'
    CSV_HEADERS = ['Genome Name',
                   'Strain',
                   'Genome Status',
                   'KEGG ID',
                   'NCBI Taxon ID',
                   'Project ID',
                   'RefSeq Project ID',
                   'Super Kingdom', 'Genus',
                   'Gram Stain', 'Shape', 'Motility', 'Pathogenic in',
                   'Genome Size', 'GC Content', 'Salinity', 'Temperature Range',
                   'Habitat', OXY_REQ,
                   'Energy Source', 'Energy Category', 'Metabolism',
                   'Sequencing Center']
    CSV_HEADER_MAPPING = dict((util.slugify(h), h) for h in CSV_HEADERS)
    CSV_HEADER_MAPPING['phylogenetic_group'] = 'Group'
    CSV_HEADER_MAPPING['phylogenetic_order'] = 'Order'
    ORG_TABLE_HEADERS = map(util.slugify, CSV_HEADERS) + ['broad_oxygen_requirement']
    ENZ_TABLE_HEADERS = ['organism', 'EC']
    
    def __init__(self, db_filename):
        self.db = database.SqliteDatabase(db_filename)
    
    def _InitTables(self):
        self.db.CreateTable('organisms', self.ORG_TABLE_HEADERS, drop_if_exists=False)
        self.db.CreateTable('organism_enzymes', self.ENZ_TABLE_HEADERS, drop_if_exists=False)
    
    def OrganismsForEC(self, ec):
        q = self.db.Execute("SELECT organism FROM organism_enzymes WHERE EC='%s'" % ec)
        for i in q:
            yield i[0]
    
    def KEGG2NCBI(self, kegg_id):
        q = self.db.Execute("SELECT ncbi_taxon_id FROM organisms WHERE kegg_id='%s'" % kegg_id)
        q = list(q)
        if not q:
            return None
        return q[0][0]
    
    def NCBI2EnergyCategory(self, ncbi_taxon):
        q = self.db.Execute("SELECT energy_category FROM organisms WHERE ncbi_taxon_id='%s'" % ncbi_taxon)
        q = list(q)
        if not q:
            return None
        return q[0][0]

    def KEGG2BroadOxygenReq(self, kegg_id):
        q = self.db.Execute("SELECT broad_oxygen_requirement FROM organisms WHERE kegg_id='%s'" % kegg_id)
        q = list(q)
        if not q:
            return None
        return q[0][0]

    def KEGG2EnergyCategory(self, kegg_id):
        q = self.db.Execute("SELECT energy_category FROM organisms WHERE kegg_id='%s'" % kegg_id)
        q = list(q)
        if not q:
            return None
        return q[0][0]
    
    def KEGG2EnergySource(self, kegg_id):
        q = self.db.Execute("SELECT energy_source FROM organisms WHERE kegg_id='%s'" % kegg_id)
        q = list(q)
        if not q:
            return None
        return q[0][0]
    
    def KEGG2Metabolism(self, kegg_id):
        q = self.db.Execute("SELECT metabolism FROM organisms WHERE kegg_id='%s'" % kegg_id)
        q = list(q)
        if not q:
            return None
        return q[0][0]
    
    @staticmethod
    def GetBroadyOxyReq(req_str):
        if not req_str:
            return None
        
        l_req = req_str.lower()        
        if 'anaer' in l_req:
            return 'anaerobe'
        if 'facult' in l_req:
            return 'facultative'
        if 'microaero' in l_req:
            return 'microaerophile'
        if 'aerob':
            return 'aerobe'
        
        raise ValueError('Couldnt Parse!')
    
    def Populate(self, filename):
        """Populates the database from files."""
        self._InitTables()
        
        f = open(filename)
        r = csv.DictReader(f)
        
        for row in r:
            insert_row = []
            for table_header in self.ORG_TABLE_HEADERS:
                if table_header not in self.CSV_HEADER_MAPPING:
                    insert_row.append(None)
                    continue
                
                csv_header = self.CSV_HEADER_MAPPING[table_header]
                val = row.get(csv_header, None)
                if val and val.strip():
                    insert_row.append(val)
                else: 
                    insert_row.append(None)
                
            oxy_req = row.get(self.OXY_REQ, None)
            broad_req = self.GetBroadyOxyReq(oxy_req)
            insert_row[-1] = broad_req
            
            self.db.Insert('organisms', insert_row)     
        f.close()
        
        k = Kegg.getInstance(loadFromAPI=False)
        enzyme_map = k.ec2enzyme_map
        for ec, enzyme in enzyme_map.iteritems():
            for org in enzyme.genes.keys():                
                self.db.Insert('organism_enzymes', [org.lower(), ec])
    
        
if __name__ == '__main__':
    db = GenomeDB('../res/genomes.sqlite')
    db.Populate('../data/genomics/IMG_tagged_genomes.csv')