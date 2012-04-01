from SOAPpy import WSDL
from toolbox.database import SqliteDatabase
import sys

class KeggGenes(object):
    
    def __init__(self):
        self.wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        self.serv = WSDL.Proxy(self.wsdl)

        self.db = SqliteDatabase('../res/channeling.sqlite', 'w')
        
        self.GENE_TABLE_NAME = 'kegg_genes'
        self.ENZYME_TABLE_NAME = 'kegg_enzymes'
        self.REACTION_TABLE_NAME = 'kegg_enzymes'
        
        self.db.CreateTable(self.GENE_TABLE_NAME, ['organism', 'gene'], drop_if_exists=False)
        self.db.CreateTable(self.ENZYME_TABLE_NAME, ['organism', 'gene', 'enzyme'], drop_if_exists=False)
        self.db.CreateTable(self.REACTION_TABLE_NAME, ['enzyme', 'reaction'], drop_if_exists=False)
    
    def GetAllGenes(self, organism='eco'):
        all_genes = []
        i = 1
        while True:
            sys.stderr.write('reading genes %d-%d ...\n' % (i, i+100))
            new_genes = self.serv.get_genes_by_organism('eco', i, 100)
            if len(new_genes) == 0:
                break
            i += 100
            all_genes += new_genes

        # clear all entries for this organism before reinserting them into the DB        
        self.db.Execute("DELETE * FROM %s WHERE organism = '%s'" % 
                        (self.GENE_TABLE_NAME, organism))
        for gene in all_genes:
            self.db.Insert(self.GENE_TABLE_NAME, [organism, gene])
        self.db.Commit()
    
    def GetAllEnzyme(self, organism='eco'):
        all_genes = []
        for row in self.db.Execute("SELECT gene FROM %s WHERE organism = '%s'" % 
                                   (self.GENE_TABLE_NAME, organism)):
            all_genes.append(str(row[0]))

        self.db.Execute("DELETE FROM %s WHERE organism = '%s'" % 
                        (self.ENZYME_TABLE_NAME, organism))
        
        for gene in all_genes:
            sys.stderr.write('reading enzymes for gene %s ...\n' % (gene))
            new_enzymes = self.serv.get_enzymes_by_gene(gene)
            for enzyme in new_enzymes:
                self.db.Insert(self.ENZYME_TABLE_NAME, [organism, gene, enzyme])
        self.db.Commit()
        
    def GetAllReactions(self):
        all_enzymes = []
        for row in self.db.Execute("SELECT distinct(enzyme) FROM %s" % 
                                   (self.ENZYME_TABLE_NAME)):
            all_enzymes.append(str(row[0]))
        
        self.db.Execute("DELETE FROM %s" % (self.REACTION_TABLE_NAME))
        
        for enzyme in all_enzymes:
            sys.stderr.write('reading reaction for enzyme %s ...\n' % (enzyme))
            new_reactions = self.serv.get_reactions_by_enzyme(enzyme)
            for reaction in new_reactions:
                self.db.Insert(self.REACTION_TABLE_NAME, [enzyme, reaction])
        self.db.Commit()

if __name__ == "__main__":
    kegg = KeggGenes()
    #kegg.GetAllGenes('eco')
    #kegg.GetAllEnzyme('eco')
    kegg.GetAllReactions()
    
