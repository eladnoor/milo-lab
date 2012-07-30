import numpy as np
import matplotlib.pyplot as plt
from toolbox.database import SqliteDatabase
import sys, csv
from pygibbs import kegg_parser
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggParseException, KeggNonCompoundException,\
    KeggReactionNotBalancedException
import itertools
from toolbox.plotting import cdf
import urllib
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_estimators import LoadAllEstimators

class KeggGenes(object):
    
    def __init__(self, html_fname):
        self.serv = None
        self.db = SqliteDatabase('channeling/channeling.sqlite', 'w')
        self.html_writer = HtmlWriter(html_fname)
        
        self.COMPOUND_TABLE_NAME = 'kegg_compounds'
        self.GENE_TABLE_NAME = 'kegg_genes'
        self.GENE_REACTION_TABLE_NAME = 'kegg_genes_to_reactions'
        self.REACTION_TABLE_NAME = 'kegg_reactions'
        self.EQUATION_TABLE_NAME = 'kegg_equations'
        self.STOICHIOMETRY_TABLE_NAME = 'kegg_stoichiometry'
        self.GIBBS_ENERGY_TABLE_NAME = 'kegg_gibbs_energies'
        self.GENE_ENERGY_TABLE_NAME = 'kegg_gene_energies'
        self.FUNCTIONAL_INTERATCTIONS_TABLE = 'parkinson_functional_interactions'
        self.GENE_PAIRS_TABLE_NAME = 'kegg_gene_pairs'
        self.COFACTOR_TABLE_NAME = 'kegg_cofactors'

    def GetAllCompounds(self):
        self.db.CreateTable(self.COMPOUND_TABLE_NAME, "compound INT, name TEXT, all_names TEXT", drop_if_exists=True)
        self.db.CreateIndex('compound_idx', self.COMPOUND_TABLE_NAME, 'compound', unique=True, drop_if_exists=True)

        f = urllib.urlopen('http://rest.kegg.jp/list/cpd/')
        for row in f.read().split('\n'):
            compound, all_names = row.split('\t')
            name = all_names.split(';')[0]
            self.db.Insert(self.COMPOUND_TABLE_NAME, [compound, name, all_names])
        self.db.Commit()
    
    def GetAllGenes(self, organism='eco'):
        self.db.CreateTable(self.GENE_TABLE_NAME, ['organism', 'gene', 'desc'], drop_if_exists=False)
        self.db.CreateIndex('gene_idx', self.GENE_TABLE_NAME, 'gene', unique=False, drop_if_exists=False)

        self.db.Execute("DELETE FROM %s WHERE organism = '%s'" % 
                        (self.GENE_TABLE_NAME, organism))

        f = urllib.urlopen('http://rest.kegg.jp/list/%s/' % organism)
        for row in f.read().split('\n'):
            gene, desc = row.split('\t')
            self.db.Insert(self.GENE_TABLE_NAME, [organism, gene, desc])
        self.db.Commit()
    
    def GetAllReactions(self, organism='eco'):
        self.db.CreateTable(self.GENE_REACTION_TABLE_NAME, ['organism', 'gene', 'reaction'], drop_if_exists=False)
        self.db.CreateIndex('reaction_gene_idx', self.GENE_REACTION_TABLE_NAME, 'gene', unique=False, drop_if_exists=False)
        self.db.CreateIndex('reaction_idx', self.GENE_REACTION_TABLE_NAME, 'reaction', unique=False, drop_if_exists=False)

        self.db.Execute("DELETE FROM %s WHERE organism = '%s'" % 
                        (self.GENE_REACTION_TABLE_NAME, organism))

        f = urllib.urlopen('http://rest.kegg.jp/link/rn/%s' % organism)
        for row in f.read().split('\n'):
            gene, reaction = row.split('\t')
            self.db.Insert(self.GENE_REACTION_TABLE_NAME, [organism, gene, reaction])
        self.db.Commit()        
                
    def GetAllEquations(self):
        self.db.CreateTable(self.EQUATION_TABLE_NAME, ['reaction', 'equation'], drop_if_exists=True)
        self.db.CreateIndex('equation_reaction_idx', self.EQUATION_TABLE_NAME, 'reaction', unique=False, drop_if_exists=True)
        self.db.CreateIndex('equation_idx', self.EQUATION_TABLE_NAME, 'equation', unique=False, drop_if_exists=True)

        all_reactions = []
        for row in self.db.Execute("SELECT distinct(reaction) FROM %s" % 
                                   (self.GENE_REACTION_TABLE_NAME)):
            all_reactions.append(str(row[0]))
        
        for reaction in all_reactions:
            f = urllib.urlopen('http://rest.kegg.jp/get/%s' % reaction)
            for equation in self._ReadReactionEntries(f.read()):
                self.db.Insert(self.EQUATION_TABLE_NAME,
                    [reaction, equation])
                sys.stderr.write('Equation for reaction %s: %s\n' % (reaction, equation))
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
        self.db.CreateTable(self.STOICHIOMETRY_TABLE_NAME, "equation TEXT, compound TEXT, coefficient REAL", drop_if_exists=True)
        self.db.CreateIndex('stoichiometry_equation_idx', self.STOICHIOMETRY_TABLE_NAME, 'equation', unique=False, drop_if_exists=True)
        self.db.CreateIndex('stoichiometry_compound_idx', self.STOICHIOMETRY_TABLE_NAME, 'compound', unique=False, drop_if_exists=True)

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
        
        for i, equation in enumerate(all_equations):
            for compound, coefficient in all_kegg_reactions[i].iteritems():
                self.db.Insert(self.STOICHIOMETRY_TABLE_NAME,
                               [equation, "cpd:C%05d" % compound, coefficient])
    
        self.db.Commit()

    def GetForamtionEnergies(self, thermo):
        self.db.CreateTable(self.GIBBS_ENERGY_TABLE_NAME, "equation TEXT, dG0 REAL, dGc REAL", drop_if_exists=True)
        self.db.CreateIndex('gibbs_equation_idx', self.GIBBS_ENERGY_TABLE_NAME, 'equation', unique=True, drop_if_exists=True)

        all_equations = set()
        for row in self.db.Execute("SELECT distinct(equation) FROM %s" % 
                                   (self.EQUATION_TABLE_NAME)):
            all_equations.add(str(row[0]))
        
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
    
    def LoadCofactors(self):
        self.db.CreateTable(self.COFACTOR_TABLE_NAME,
                            'compound TEXT, name TEXT, c_min REAL, c_max REAL, ref TEXT',
                            drop_if_exists=True)
        self.db.CreateIndex('cofactor_idx', self.COFACTOR_TABLE_NAME,
                            'compound', unique=True, drop_if_exists=True)

        csv_reader = csv.DictReader(open('channeling/cofactors.csv', 'r'))
        for rowdict in csv_reader:
            self.db.Insert(self.COFACTOR_TABLE_NAME,
                           ["cpd:C%05d" % int(rowdict['cid']), rowdict['name'],
                            float(rowdict['c_min'] or np.nan), float(rowdict['c_max'] or np.nan),
                            rowdict['ref']])
        self.db.Commit()
        
    def CreateGeneEnergyTable(self):
        self.db.CreateTable(self.GENE_ENERGY_TABLE_NAME,
                            "gene TEXT, reaction TEXT, dGc REAL, compound INT, coefficient REAL",
                            drop_if_exists=True)
        self.db.CreateIndex('gene_energy_compound_idx',
                            self.GENE_ENERGY_TABLE_NAME, 'compound', unique=False)
        self.db.CreateIndex('gene_energy_gene_idx',
                            self.GENE_ENERGY_TABLE_NAME, 'gene', unique=False)

        query = """
            INSERT INTO %s (gene, reaction, dGc, compound, coefficient)
                SELECT  gen.gene, rxn.reaction, eng.dGc, sto.compound, sto.coefficient
                FROM    kegg_genes gen, kegg_genes_to_reactions rxn,
                        kegg_equations eqn, kegg_gibbs_energies eng,
                        kegg_stoichiometry sto
                WHERE   gen.organism = 'eco'
                AND     gen.gene = rxn.gene
                AND     rxn.reaction = eqn.reaction
                AND     eqn.equation = eng.equation
                AND     eng.dG0 IS NOT NULL
                AND     eqn.equation = sto.equation
        """ % self.GENE_ENERGY_TABLE_NAME
        self.db.Execute(query)

        query = """
            INSERT INTO %s (gene, reaction, dGc, compound, coefficient)
                SELECT  gen.gene, rxn.reaction, -eng.dGc, sto.compound, -sto.coefficient
                FROM    kegg_genes gen, kegg_genes_to_reactions rxn,
                        kegg_equations eqn, kegg_gibbs_energies eng,
                        kegg_stoichiometry sto
                WHERE   gen.organism = 'eco'
                AND     gen.gene = rxn.gene
                AND     rxn.reaction = eqn.reaction
                AND     eqn.equation = eng.equation
                AND     eng.dG0 IS NOT NULL
                AND     eqn.equation = sto.equation
        """ % self.GENE_ENERGY_TABLE_NAME
        self.db.Execute(query)
        self.db.Commit()
        
    def CreateGenePairsTable(self):
        self.db.CreateTable(self.GENE_PAIRS_TABLE_NAME,
                            "gene1 TEXT, gene2 TEXT, reaction1 TEXT, reaction2 TEXT, "
                            "compound TEXT, coeff1 REAL, coeff2 REAL, dGc1 REAL, "
                            "dGc2 REAL, score REAL",
                            drop_if_exists=True)
        self.db.CreateIndex('gene_pairs_gene_idx',
                            self.GENE_PAIRS_TABLE_NAME,
                            'gene1, gene2', unique=False)
        query = """
            INSERT INTO %s (gene1, gene2, reaction1, reaction2, compound, coeff1, coeff2, dGc1, dGc2, score)
                SELECT p.*, pfi.score FROM
                (
                    SELECT  kge1.gene gene1,
                            kge2.gene gene2, 
                            kge1.reaction reaction1,
                            kge2.reaction reaction2,
                            kge1.compound compound, 
                            kge1.coefficient coeff1,
                            kge2.coefficient coeff2,
                            cast(kge1.dGc as real) dGc1, 
                            cast(kge2.dGc as real) dGc2
                    FROM    kegg_gene_energies kge1, kegg_gene_energies kge2
                    WHERE   kge1.compound = kge2.compound
                    AND     kge1.compound NOT IN (SELECT compound FROM %s)
                    AND     kge1.gene != kge2.gene
                    AND     kge1.reaction != kge2.reaction
                    AND     kge1.coefficient > 0
                    AND     kge2.coefficient < 0
                ) p
                LEFT OUTER JOIN %s pfi
                ON      (pfi.gene1 = p.gene1 AND pfi.gene2 = p.gene2
                         OR
                         pfi.gene1 = p.gene2 AND pfi.gene2 = p.gene1)
        """ % (self.GENE_PAIRS_TABLE_NAME, self.COFACTOR_TABLE_NAME, self.FUNCTIONAL_INTERATCTIONS_TABLE)
        self.db.Execute(query)
        self.db.Commit()
        
    def Correlate(self, dGc1_lower, dGc2_upper, reverse=False):
        if reverse:
            cond = "kgp.dGc1 < %d AND kgp.dGc2 > %d" % (dGc1_lower, dGc2_upper)
        else:
            cond = "kgp.dGc1 > %d AND kgp.dGc2 < %d" % (dGc1_lower, dGc2_upper)
        
        query = """
            SELECT  kgp.gene1, kgp.gene2, sum(%s) nqual, count(*) ntot, max(score)
            FROM %s kgp
            GROUP BY kgp.gene1, kgp.gene2
        """ % (cond, self.GENE_PAIRS_TABLE_NAME)
        
        counters = np.zeros((2, 2))
        
        for row in self.db.Execute(query):
            _gene1, _gene2, nqual, _ntot, score = row
            i = int(score is not None) # is there an PP-interaction
            j = int(nqual > 0) # is this a qualifying pair (thermodynamically)
            counters[i, j] += 1.0

        _inter0 = np.sum(counters[0, :])
        inter1 = np.sum(counters[1, :])
        qual0 = np.sum(counters[:, 0])
        qual1 = np.sum(counters[:, 1])
        total = np.sum(counters.flat)
        
        print "-" * 50
        if reverse:
            print "Checking criterion: first < %d and second > %d" % (dGc1_lower, dGc2_upper)
        else:
            print "Checking criterion: first > %d and second < %d" % (dGc1_lower, dGc2_upper)
        print "Total no. of pairs = %d" % total
        
        print "interaction rate among all pairs (%d out of %d) = %.2f%%" % (inter1, total, 100*(inter1 / total))
        print "qualification rate among all pairs (%d out of %d) = %.2f%%" % (qual1, total, 100*(qual1 / total))
        print "interactions between unqualifying pairs (%d out of %d) = %.2f%%" % (counters[1,0], qual0, 100*(counters[1,0] / qual0))
        print "interactions between qualifying pairs (%d out of %d) = %.2f%%" % (counters[1,1], qual1, 100*(counters[1,1] / qual1))
        
        return counters[1,0] / qual0, counters[1,1] / qual1


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

    def PlotScatter(self):
        query = """
                SELECT  p.g1, p.g2, pfi.score
                FROM (
                      SELECT  kgp.gene1 gene1, kgp.gene2 gene2, cast(kgp.dGc1 as real) g1, cast(kgp.dGc2 as real) g2
                      FROM    %s kgp
                     ) p
                LEFT OUTER JOIN %s pfi
                ON      (pfi.gene1 = p.gene1 AND pfi.gene2 = p.gene2
                         OR
                         pfi.gene1 = p.gene2 AND pfi.gene2 = p.gene1)
            """ % (self.GENE_PAIRS_TABLE_NAME, self.FUNCTIONAL_INTERATCTIONS_TABLE)

        data = []
        for row in self.db.Execute(query):
            g1, g2, score = row
            data.append([float(g1), float(g2), float(score or 0)])
        data = np.matrix(data)

        ind1 = list(np.where(data[:, 2] > 0)[0].flat)
        ind2 = list(np.where(data[:, 2] == 0)[0].flat)
        fig = plt.figure(figsize=(6,6), dpi=90)    
        plt.plot(data[ind2, 0], data[ind2, 1], 'r.', markersize=5, figure=fig)
        plt.plot(data[ind1, 0], data[ind1, 1], 'g.', markersize=5, figure=fig)
        plt.show()
        
    def PlotCDF(self):
        special_pairs = {('eco:b3236', 'eco:b0720'):"mdh:gltA", # malate dehydrogenase -> oxaloacetate -> citrate synthase
                         ('eco:b1263', 'eco:b1264'):"trpD:trpE"} # trpD -> chorismate -> trpE (two components of anthraline synthase)
        
        query = """
                SELECT gene1, gene2, min(dGc2 - dGc1), max(score)
                FROM %s
                WHERE dGc1 + dGc2 < 0
                AND dGc1 > 10
                GROUP BY gene1, gene2
                """ % (self.GENE_PAIRS_TABLE_NAME)

        data = []
        markers = []
        for row in self.db.Execute(query):
            gene1, gene2, ddG, score = row
            if (gene1, gene2) in special_pairs:
                markers.append((special_pairs[(gene1, gene2)], ddG))
            data.append([ddG, float(score or 0)])
        data = np.matrix(data)

        ind1 = list(np.where(data[:, 1] > 0)[0].flat)
        ind2 = list(np.where(data[:, 1] == 0)[0].flat)
    
        fig = plt.figure(figsize=(6,6), dpi=90)    
        cdf((data[ind2, 0]).flat, label="non-interacting (N = %d)" % len(ind2), style='r', figure=fig)
        cdf((data[ind1, 0]).flat, label="interacting (N = %d)" % len(ind1), style='g', figure=fig)
        for label, ddG in markers:
            plt.plot([ddG, ddG], [0, 1], 'b--', figure=fig)
            plt.text(ddG, 0.1, label)
        plt.xlim(-500, 500)
        plt.xlabel(r"$\Delta G'^c$ (2nd) - $\Delta G'^c$ (1st) [kJ/mol]")
        plt.ylabel(r"Cumulative Distribution Function")
        plt.legend(loc="upper left")

        self.html_writer.embed_matplotlib_figure(fig, width=400, height=400, name='channeling_cdf')

    def PrintPairs(self):
        query = """
                SELECT g.gene1, g.gene2, c.name, g.reaction1, g.reaction2, 
                       cast(g.dG1 as int), cast(g.dG2 as int), cast(g.ddG as int),
                       kg1.desc, kg2.desc, g.score FROM
                (SELECT gene1, gene2, reaction1, reaction2, compound, max(dGc1) dG1, min(dGc2) dG2, min(dGc2 - dGc1) ddG, max(score) score
                FROM kegg_gene_pairs
                WHERE dGc1 + dGc2 < 0
                AND   dGc1 > 10
                GROUP BY gene1, gene2, compound
                ORDER BY ddG) g, kegg_genes kg1, kegg_genes kg2, kegg_compounds c
                WHERE g.gene1 = kg1.gene AND g.gene2 = kg2.gene AND c.compound = g.compound
                """
        
        self.html_writer.write('<font size="1">\n')
        self.db.Query2HTML(self.html_writer, query,
                           ['Gene 1', 'Gene 2', 'Common Compound',
                            'Reaction 1', 'Reaction 2',
                            'dGc1', 'dGc2', 'dG2-dG1', 'Desc 1', 'Desc 2',
                            'Score'])
        self.html_writer.write('</font>\n')

if __name__ == "__main__":
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 2
    
    kegg_gene = KeggGenes('../res/channeling.html')

    if False:
        kegg_gene.LoadCofactors()
        kegg_gene.LoadFunctionalInteractions()
        kegg_gene.GetAllCompounds()
        kegg_gene.GetAllGenes('eco')
        kegg_gene.GetAllReactions('eco')
        kegg_gene.GetAllEquations()
        kegg_gene.GetStoichiometries()
        
        estimators = LoadAllEstimators()
        kegg_gene.GetForamtionEnergies(estimators['UGC'])
        
        kegg_gene.LoadCofactors()
        kegg_gene.CreateGeneEnergyTable()
        kegg_gene.CreateGenePairsTable()
    
    #kegg_gene.Correlate(-40, -40, reverse=False)
    #kegg_gene.Correlate(-40, -30, reverse=False)
    kegg_gene.PlotCDF()
    kegg_gene.PrintPairs()
    sys.exit(0)
    
    csv_writer = csv.writer(open('../res/channeling.csv', 'w'))
    csv_writer.writerow(['dGc1 <> x', 'dGc2 <> x', 'P(unqualify)', 'P(qualify)'])
    energies1 = np.arange(-40, 41, 10)
    energies2 = np.arange(-40, 41, 10)
    
    for e1, e2 in itertools.product(energies1, energies2):
        p0, p1 = kegg_gene.Correlate(e1, e2, reverse=False)
        csv_writer.writerow([" > %g" % e1, " < %g" % e2, p0, p1])

    for e1, e2 in itertools.product(energies1, energies2):
        p0, p1 = kegg_gene.Correlate(e1, e2, reverse=True)
        csv_writer.writerow([" < %g" % e1, " > %g" % e2, p0, p1])
    