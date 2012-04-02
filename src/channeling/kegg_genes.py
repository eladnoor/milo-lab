from SOAPpy import WSDL
from toolbox.database import SqliteDatabase
import sys
from pygibbs import kegg_parser
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggParseException, KeggNonCompoundException

class KeggGenes(object):
    
    def __init__(self):
        self.wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        self.serv = WSDL.Proxy(self.wsdl)

        self.db = SqliteDatabase('channeling/channeling.sqlite', 'w')
        
        self.GENE_TABLE_NAME = 'kegg_genes'
        self.ENZYME_TABLE_NAME = 'kegg_enzymes'
        self.REACTION_TABLE_NAME = 'kegg_reactions'
        self.EQUATION_TABLE_NAME = 'kegg_equations'
        self.GIBBS_ENERGY_TABLE_NAME = 'kegg_gibbs_energies'
        
        self.db.CreateTable(self.GENE_TABLE_NAME, ['organism', 'gene'], drop_if_exists=False)
        self.db.CreateTable(self.ENZYME_TABLE_NAME, ['organism', 'gene', 'enzyme'], drop_if_exists=False)
        self.db.CreateTable(self.REACTION_TABLE_NAME, ['enzyme', 'reaction'], drop_if_exists=False)
        self.db.CreateTable(self.EQUATION_TABLE_NAME, ['reaction', 'equation'], drop_if_exists=False)
        self.db.CreateTable(self.GIBBS_ENERGY_TABLE_NAME, ['equation', 'dG0'], drop_if_exists=False)
    
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
            s = self.serv.bget(reaction)
            for equation in self._ReadReactionEntries(s):
                sys.stderr.write(equation + "\n")
                self.db.Insert(self.EQUATION_TABLE_NAME,
                    [reaction, equation])
        self.db.Commit()

    def _ReadReactionEntries(self, s):
        equation_list = []
        entry2fields_map = kegg_parser.ParsedKeggFile.FromKeggAPI(s)
        for key in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[key]
            if "EQUATION" in field_map:
                equation_list.append(field_map["EQUATION"])
        
        return equation_list
    
    def GetForamtionEnergies(self, thermo):
        from pygibbs.kegg import Kegg
        kegg = Kegg.getInstance()
        all_kegg_cids = set(kegg.get_all_cids())
    
        all_kegg_reactions = []
        all_equations = []
        for row in self.db.Execute("SELECT distinct(equation) FROM %s" % 
                                   (self.EQUATION_TABLE_NAME)):
            try:
                r = Reaction.FromFormula(str(row[0]))
                if not r.get_cids().issubset(all_kegg_cids):
                    continue
                all_equations.append(str(row[0]))
                all_kegg_reactions.append(r)
            except (KeggParseException, KeggNonCompoundException):
                pass
        
        self.db.Execute("DELETE FROM %s" % (self.GIBBS_ENERGY_TABLE_NAME))
        
        dG0 = thermo.GetTransfromedKeggReactionEnergies(all_kegg_reactions)
        for i, equation in enumerate(all_equations):
            sys.stderr.write('estimating energy for %s ...\n' % str(all_kegg_reactions[i]))
            self.db.Insert(self.GIBBS_ENERGY_TABLE_NAME,
                           [equation, dG0[0, i]])
    
        self.db.Commit()

if __name__ == "__main__":
    kegg_gene = KeggGenes()
    #kegg_gene.GetAllGenes('eco')
    #kegg_gene.GetAllEnzyme('eco')
    #kegg_gene.GetAllReactions()
    #kegg_gene.GetAllEquations()
    
    from pygibbs.thermodynamic_estimators import LoadAllEstimators
    estimators = LoadAllEstimators()
    kegg_gene.GetForamtionEnergies(estimators['UGC'])
    
    
