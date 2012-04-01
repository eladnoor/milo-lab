from SOAPpy import WSDL
from toolbox.database import SqliteDatabase
from pygibbs.kegg import Kegg
import sys

class KeggGenes(object):
    
    def __init__(self):
        self.kegg = Kegg.getInstance()
    
        self.wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        self.serv = WSDL.Proxy(self.wsdl)

        self.db = SqliteDatabase('channeling/channeling.sqlite', 'w')
        
        self.GENE_TABLE_NAME = 'kegg_genes'
        self.ENZYME_TABLE_NAME = 'kegg_enzymes'
        self.REACTION_TABLE_NAME = 'kegg_reactions'
        self.EQUATION_TABLE_NAME = 'kegg_equations'
        
        self.db.CreateTable(self.GENE_TABLE_NAME, ['organism', 'gene'], drop_if_exists=False)
        self.db.CreateTable(self.ENZYME_TABLE_NAME, ['organism', 'gene', 'enzyme'], drop_if_exists=False)
        self.db.CreateTable(self.REACTION_TABLE_NAME, ['enzyme', 'reaction'], drop_if_exists=False)
        self.db.CreateTable(self.EQUATION_TABLE_NAME, ['reaction', 'equation', 'dG0'], drop_if_exists=False)
    
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
            sys.stderr.write('reading reactions for enzyme %s ...\n' % (enzyme))
            new_reactions = self.serv.get_reactions_by_enzyme(enzyme)
            for reaction in new_reactions:
                self.db.Insert(self.REACTION_TABLE_NAME, [enzyme, reaction])
        self.db.Commit()

    def GetAllEquations(self):
        all_reactions = []
        for row in self.db.Execute("SELECT distinct(reaction) FROM %s" % 
                                   (self.REACTION_TABLE_NAME)):
            all_reactions.append(str(row[0]))
        
        self.db.Execute("DELETE FROM %s" % (self.EQUATION_TABLE_NAME))
        
        for reaction in all_reactions:
            sys.stderr.write('reading data for reaction %s ...\n' % (reaction))
            rid = int(reaction[4:])
            try:
                kegg_reaction = self.kegg.rid2reaction(rid)
                equation = kegg_reaction.equation
                self.db.Insert(self.EQUATION_TABLE_NAME,
                    [reaction, equation, None])
            except KeyError:
                sys.stderr.write("Cannot find this RID in local KEGG copy")
        self.db.Commit()
        

if __name__ == "__main__":
    kegg_gene = KeggGenes()
    #kegg_gene.GetAllGenes('eco')
    #kegg_gene.GetAllEnzyme('eco')
    #kegg_gene.GetAllReactions()
    kegg_gene.GetAllEquations()
    
    
