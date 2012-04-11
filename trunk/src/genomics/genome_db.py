import csv

from toolbox import database
from toolbox import util

class GenomeDB(object):
    
    CSV_HEADERS = ['Organism Name',
                   'Strain',
                   'Genome Status',
                   'KEGG ID',
                   'NCBI Taxon ID',
                   'Project ID',
                   'RefSeq Project ID',
                   'Super Kingdom', 'Genus',
                   'Gram Stain', 'Shape', 'Motility', 'Pathogenic in',
                   'Genome Size', 'GC Content', 'Salinity', 'Temperature Range',
                   'Habitat', 'Oxygen Requirement IMG', 'Energy Source',
                   'Energy Category', 'Metabolism', 'Sequencing Center']
    CSV_HEADER_MAPPING = dict((util.slugify(h), h) for h in CSV_HEADERS)
    CSV_HEADER_MAPPING['phylogenetic_group'] = 'Group'
    CSV_HEADER_MAPPING['phylogenetic_order'] = 'Order'
    TABLE_HEADERS = map(util.slugify, CSV_HEADERS) + ['broad_oxygen_requirement']
    
    def __init__(self, db_filename):
        self.db = database.SqliteDatabase(db_filename)
        self._InitTables()
    
    def _InitTables(self):
        self.db.CreateTable('organisms', self.TABLE_HEADERS, drop_if_exists=False)
    
    def PopulateFromFile(self, filename):
        f = open(filename)
        r = csv.DictReader(f)
        
        for row in r:
            insert_row = []
            for table_header in self.TABLE_HEADERS:
                if table_header not in self.CSV_HEADER_MAPPING:
                    insert_row.append(None)
                    continue
                
                csv_header = self.CSV_HEADER_MAPPING[table_header]
                insert_row.append(row.get(csv_header, None))
                
            self.db.Insert('organisms', insert_row)     
        f.close()
    
    
        
if __name__ == '__main__':
    db = GenomeDB('../res/genomes.sqlite')
    db.PopulateFromFile('../data/genomics/merged_tagged_sequenced_genomes.csv')