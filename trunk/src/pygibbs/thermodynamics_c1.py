from pygibbs.thermodynamics import ReactionThermodynamics
import csv
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg import Kegg

def GetC1Thermodynamics(reaction_fname='../data/thermodynamics/c1_reaction_thermodynamics.csv'):
    kegg = Kegg.getInstance()
    for row in csv.DictReader(open(reaction_fname, 'r')):
        r = Reaction.FromFormula(row['formula'])
        r.SetNames(row['enzyme'])
    