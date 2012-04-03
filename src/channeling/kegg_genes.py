from SOAPpy import WSDL
from toolbox.database import SqliteDatabase
import sys, csv
from pygibbs import kegg_parser
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggParseException, KeggNonCompoundException,\
    KeggReactionNotBalancedException
import itertools

class KeggGenes(object):
    
    def __init__(self):
        self.db = SqliteDatabase('channeling/channeling.sqlite', 'w')
        
        self.KEGG_WSDL = 'http://soap.genome.jp/KEGG.wsdl'
        self.GENE_TABLE_NAME = 'kegg_genes'
        self.ENZYME_TABLE_NAME = 'kegg_enzymes'
        self.REACTION_TABLE_NAME = 'kegg_reactions'
        self.EQUATION_TABLE_NAME = 'kegg_equations'
        self.STOICHIOMETRY_TABLE_NAME = 'kegg_stoichiometry'
        self.GIBBS_ENERGY_TABLE_NAME = 'kegg_gibbs_energies'
        self.GENE_ENERGY_TABLE_NAME = 'kegg_gene_energies'
        self.FUNCTIONAL_INTERATCTIONS_TABLE = 'parkinson_functional_interactions'
        self.GENE_PAIRS_TABLE_NAME = 'kegg_gene_pairs'
        
        self.db.CreateTable(self.GENE_TABLE_NAME, ['organism', 'gene'], drop_if_exists=False)
        self.db.CreateIndex('gene_idx', self.GENE_TABLE_NAME, 'gene', unique=False)
        
        self.db.CreateTable(self.ENZYME_TABLE_NAME, ['organism', 'gene', 'enzyme'], drop_if_exists=False)
        self.db.CreateIndex('enzyme_gene_idx', self.ENZYME_TABLE_NAME, 'gene', unique=False)
        self.db.CreateIndex('enzyme_idx', self.ENZYME_TABLE_NAME, 'enzyme', unique=False)

        self.db.CreateTable(self.REACTION_TABLE_NAME, ['enzyme', 'reaction'], drop_if_exists=False)
        self.db.CreateIndex('reaction_enzyme_idx', self.REACTION_TABLE_NAME, 'enzyme', unique=False)
        self.db.CreateIndex('reaction_idx', self.REACTION_TABLE_NAME, 'reaction', unique=False)
        
        self.db.CreateTable(self.EQUATION_TABLE_NAME, ['reaction', 'equation'], drop_if_exists=False)
        self.db.CreateIndex('equation_reaction_idx', self.EQUATION_TABLE_NAME, 'reaction', unique=False)
        self.db.CreateIndex('equation_idx', self.EQUATION_TABLE_NAME, 'equation', unique=False)
        
        self.db.CreateTable(self.STOICHIOMETRY_TABLE_NAME, ['equation', 'compound', 'coefficient'], drop_if_exists=False)
        self.db.CreateIndex('stoichiometry_equation_idx', self.STOICHIOMETRY_TABLE_NAME, 'equation', unique=False)
        self.db.CreateIndex('stoichiometry_compound_idx', self.STOICHIOMETRY_TABLE_NAME, 'compound', unique=False)

        self.db.CreateTable(self.GIBBS_ENERGY_TABLE_NAME, ['equation', 'dG0', 'dGc'], drop_if_exists=False)
        self.db.CreateIndex('gibbs_equation_idx', self.GIBBS_ENERGY_TABLE_NAME, 'equation', unique=True)
        
    def get_kegg_serv(self):
        if self.serv is None:
            self.serv = WSDL.Proxy(self.KEGG_WSDL)
        return self.serv
            
    def GetAllGenes(self, organism='eco'):
        all_genes = []
        i = 1
        while True:
            sys.stderr.write('reading genes %d-%d ...\n' % (i, i+100))
            new_genes = self.get_kegg_serv().get_genes_by_organism('eco', i, 100)
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
            new_enzymes = self.get_kegg_serv().get_enzymes_by_gene(gene)
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
            new_reactions = self.get_kegg_serv().get_reactions_by_enzyme(enzyme)
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
            s = self.get_kegg_serv().bget(reaction)
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
    
    def GetStoichiometries(self):
        all_kegg_reactions = []
        all_equations = []
        for row in self.db.Execute("SELECT distinct(equation) FROM %s" % 
                                   (self.EQUATION_TABLE_NAME)):
            try:
                r = Reaction.FromFormula(str(row[0]))
                all_equations.append(str(row[0]))
                all_kegg_reactions.append(r)
            except (KeggParseException, KeggNonCompoundException):
                pass
        
        self.db.Execute("DELETE FROM %s" % (self.STOICHIOMETRY_TABLE_NAME))
        
        for i, equation in enumerate(all_equations):
            for compound, coefficient in all_kegg_reactions[i].iteritems():
                self.db.Insert(self.STOICHIOMETRY_TABLE_NAME,
                               [equation, compound, coefficient])
    
        self.db.Commit()

    def GetForamtionEnergies(self, thermo):
        all_equations = set()
        for row in self.db.Execute("SELECT distinct(equation) FROM %s" % 
                                   (self.EQUATION_TABLE_NAME)):
            all_equations.add(str(row[0]))
        
        self.db.Execute("DELETE FROM %s" % (self.GIBBS_ENERGY_TABLE_NAME))
        from pygibbs.kegg import Kegg
        kegg = Kegg.getInstance()
        all_kegg_cids = set(kegg.get_all_cids())
        for equation in all_equations:
            try:
                rxn = Reaction.FromFormula(equation)
                if not rxn.get_cids().issubset(all_kegg_cids):
                    raise KeggNonCompoundException
                rxn.Balance(balance_water=True, exception_if_unknown=True)
                dG0 = thermo.GetTransfromedKeggReactionEnergies([rxn], conc=1)[0, 0]
                dGc = thermo.GetTransfromedKeggReactionEnergies([rxn], conc=1e-3)[0, 0]
                self.db.Insert(self.GIBBS_ENERGY_TABLE_NAME, [equation, dG0, dGc])
                
            except (KeggParseException, KeggNonCompoundException, KeggReactionNotBalancedException):
                self.db.Insert(self.GIBBS_ENERGY_TABLE_NAME, [equation, None, None])
    
        self.db.Commit()
        
    def CreateGeneEnergyTable(self):
        self.db.CreateTable(self.GENE_ENERGY_TABLE_NAME,
                            ['gene', 'enzyme', 'dGc', 'compound', 'coefficient'],
                            drop_if_exists=True)
        self.db.CreateIndex('gene_energy_compound_idx',
                            self.GENE_ENERGY_TABLE_NAME, 'compound', unique=False)
        self.db.CreateIndex('gene_energy_gene_idx',
                            self.GENE_ENERGY_TABLE_NAME, 'gene', unique=False)

        query = """
            INSERT INTO kegg_gene_energies (gene, enzyme, dGc, compound, coefficient)
                SELECT  gen.gene, enz.enzyme, eng.dGc, sto.compound, sto.coefficient
                FROM    kegg_genes gen, kegg_enzymes enz, kegg_reactions rxn,
                        kegg_equations eqn, kegg_gibbs_energies eng,
                        kegg_stoichiometry sto
                WHERE   gen.organism = 'eco'
                AND     gen.gene = enz.gene
                AND     enz.enzyme = rxn.enzyme
                AND     rxn.reaction = eqn.reaction
                AND     eqn.equation = eng.equation
                AND     eng.dG0 IS NOT NULL
                AND     eqn.equation = sto.equation
        """
        self.db.Execute(query)

        query = """
            INSERT INTO kegg_gene_energies (gene, enzyme, dGc, compound, coefficient)
                SELECT  gen.gene, enz.enzyme, -eng.dGc, sto.compound, -sto.coefficient
                FROM    kegg_genes gen, kegg_enzymes enz, kegg_reactions rxn,
                        kegg_equations eqn, kegg_gibbs_energies eng,
                        kegg_stoichiometry sto
                WHERE   gen.organism = 'eco'
                AND     gen.gene = enz.gene
                AND     enz.enzyme = rxn.enzyme
                AND     rxn.reaction = eqn.reaction
                AND     eqn.equation = eng.equation
                AND     eng.dG0 IS NOT NULL
                AND     eqn.equation = sto.equation
        """
        self.db.Execute(query)
        self.db.Commit()
        
    def CreateGenePairsTable(self):
        self.db.CreateTable(self.GENE_PAIRS_TABLE_NAME,
                            ['gene1', 'gene2', 'enzyme1', 'enzyme2', 'compound',
                             'coeff1', 'coeff2', 'dGc1', 'dGc2'],
                            drop_if_exists=True)
        self.db.CreateIndex('gene_pairs_gene_idx',
                            self.GENE_PAIRS_TABLE_NAME,
                            'gene1, gene2', unique=False)
        query = """
            INSERT INTO kegg_gene_pairs (gene1, gene2, enzyme1, enzyme2, compound, coeff1, coeff2, dGc1, dGc2)
                SELECT  kge1.gene, kge2.gene, 
                        kge1.enzyme, kge2.enzyme,
                        kge1.compound, 
                        kge1.coefficient, kge2.coefficient,
                        kge1.dGc, kge2.dGc
                FROM    kegg_gene_energies kge1, kegg_gene_energies kge2
                WHERE   kge1.compound = kge2.compound
                AND     kge1.compound IS NOT 80
                AND     kge1.compound IS NOT 1
                AND     kge1.gene != kge2.gene
                AND     kge1.enzyme != kge2.enzyme
                AND     kge1.coefficient > 0
                AND     kge2.coefficient < 0
        """
        self.db.Execute(query)
        self.db.Commit()
        
    def Correlate(self, dGc1_lower, dGc2_upper):
        queries = []
        
        queries.append("""
            SELECT  count(*), (pfi.score IS NOT NULL) s, (kgp.dGc1 > %d AND kgp.dGc2 < %d) e
            FROM    kegg_gene_pairs kgp
            LEFT OUTER JOIN parkinson_functional_interactions pfi
            ON      (pfi.gene1 = kgp.gene1 AND pfi.gene2 = kgp.gene2
                     OR
                     pfi.gene1 = kgp.gene2 AND pfi.gene2 = kgp.gene1)
            GROUP BY s, e
        """ % (dGc1_lower, dGc2_upper))
        
        res = [row for row in self.db.Execute(queries[0])]
        
        inter0_energy0 = float(res[0][0])
        inter0_energy1 = float(res[1][0])
        inter1_energy0 = float(res[2][0])
        inter1_energy1 = float(res[3][0])
        inter0 = inter0_energy0 + inter0_energy1
        inter1 = inter1_energy0 + inter1_energy1
        energy0 = inter0_energy0 + inter1_energy0
        energy1 = inter0_energy1 + inter1_energy1
        total = inter0 + inter1
        
        print "Total no. of pairs = %d" % total
        print "Total no. of interacting pairs = %d" % inter1
        print "Total no. of qualifying pairs = %d" % energy1
        
        print "interactions among all pairs (%d out of %d) = %.2f%%" % (inter1, total, 100.0 * (inter1 / total))
        print "interactions between unqualified pairs (%d out of %d) = %.2f%%" % (inter1_energy0, energy0, 100.0 * (inter1_energy0 / energy0))
        print "interactions between qualified pairs (%d out of %d) = %.2f%%" % (inter1_energy1, energy1, 100.0 * (inter1_energy1 / energy1))
        
        return inter1_energy0 / energy0, inter1_energy1 / energy1, energy1

    def Correlate2(self, dGc1_lower, dGc2_upper):
        queries = []
        
        queries.append("""
            SELECT  p.*, pfi.score
            FROM (
                  SELECT  kgp.gene1 gene1, kgp.gene2 gene2, sum(kgp.dGc1 > %d AND kgp.dGc2 < %d) nqual, count(*) ntot
                  FROM    kegg_gene_pairs kgp
                  GROUP BY kgp.gene1, kgp.gene2
                 ) p
            LEFT OUTER JOIN parkinson_functional_interactions pfi
            ON      (pfi.gene1 = p.gene1 AND pfi.gene2 = p.gene2
                     OR
                     pfi.gene1 = p.gene2 AND pfi.gene2 = p.gene1)
        """ % (dGc1_lower, dGc2_upper))
        
        inter0_energy0 = 0.0
        inter0_energy1 = 0.0
        inter1_energy0 = 0.0
        inter1_energy1 = 0.0

        for row in self.db.Execute(queries[0]):
            gene1, gene2, nqual, ntot, score = row
            if nqual == 0 and score is None:
                inter0_energy0 += 1.0
            elif nqual > 0 and score is None:
                inter0_energy1 += 1.0
            elif nqual == 0 and score is not None:
                inter1_energy0 += 1.0
            elif nqual > 0 and score is not None:
                inter1_energy1 += 1.0

        inter0 = inter0_energy0 + inter0_energy1
        inter1 = inter1_energy0 + inter1_energy1
        energy0 = inter0_energy0 + inter1_energy0
        energy1 = inter0_energy1 + inter1_energy1
        total = inter0 + inter1
        
        print "Total no. of pairs = %d" % total
        print "Total no. of interacting pairs = %d" % inter1
        print "Total no. of qualifying pairs = %d" % energy1
        
        print "interactions among all pairs (%d out of %d) = %.2f%%" % (inter1, total, 100*(inter1 / total))
        print "interactions between unqualified pairs (%d out of %d) = %.2f%%" % (inter1_energy0, energy0, 100*(inter1_energy0 / energy0))
        print "interactions between qualified pairs (%d out of %d) = %.2f%%" % (inter1_energy1, energy1, 100*(inter1_energy1 / energy1))
        
        return inter1_energy0 / energy0, inter1_energy1 / energy1, energy1


    def LoadFunctionalInteractions(self,
            fname='../data/proteomics/coli/functional_interactions.txt'):

        self.db.CreateTable(self.FUNCTIONAL_INTERATCTIONS_TABLE,
                            ['gene1', 'gene2', 'score'],
                            drop_if_exists=True)
        self.db.CreateIndex('interaction_gene_idx',
                            self.FUNCTIONAL_INTERATCTIONS_TABLE,
                            'gene1, gene2', unique=False)
        
        tsv = csv.reader(open(fname, 'r'), delimiter='\t')
        for row in tsv:
            if row[0][0] == '#':
                continue
            gene1 = 'eco:' + row[0].lower()
            gene2 = 'eco:' + row[1].lower()
            score = float(row[2])
            self.db.Insert(self.FUNCTIONAL_INTERATCTIONS_TABLE,
                           [gene1, gene2, score])
        
        self.db.Commit()

if __name__ == "__main__":
    kegg_gene = KeggGenes()
    #kegg_gene.LoadFunctionalInteractions()
    #kegg_gene.GetAllGenes('eco')
    #kegg_gene.GetAllEnzyme('eco')
    #kegg_gene.GetAllReactions()
    #kegg_gene.GetAllEquations()
    #kegg_gene.GetStoichiometries()
    
    #from pygibbs.thermodynamic_estimators import LoadAllEstimators
    #estimators = LoadAllEstimators()
    #kegg_gene.GetForamtionEnergies(estimators['UGC'])
    
    #kegg_gene.CreateGeneEnergyTable()
    #kegg_gene.CreateGenePairsTable()
    kegg_gene.Correlate2(0, 0)
    sys.exit(0)
    
    csv_writer = csv.writer(open('../res/channeling.csv', 'w'))
    csv_writer.writerow(['dGc1 > x', 'dGc2 < y', '# qualify', 'P(unqualify)', 'P(qualify)'])
    energies1 = [0, 10, 20, 30, 40]
    energies2 = [0, -10, -20, -30, -40]
    for e1, e2 in itertools.product(energies1, energies2):
        p0, p1, n_qual = kegg_gene.Correlate(e1, e2)
        csv_writer.writerow([e1, e2, n_qual, p0, p1])
    