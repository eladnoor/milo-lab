from pygibbs.thermodynamics import ReactionThermodynamics,\
    PsuedoisomerTableThermodynamics, BinaryThermodynamics
import csv
from pygibbs.kegg_reaction import Reaction
from toolbox.html_writer import HtmlWriter
from toolbox.database import SqliteDatabase

def GetC1Thermodynamics(html_writer, reaction_fname='../data/thermodynamics/c1_reaction_thermodynamics.csv'):
    html_writer.write("<h1>C1 thermodynamics</h1>\n")
    
    dict_list = []
    db_public = SqliteDatabase('../data/public_data.sqlite')
    alberty = PsuedoisomerTableThermodynamics.FromDatabase(\
                        db_public, 'alberty_pseudoisomers', name='alberty')
    alberty.AddPseudoisomer(101, nH=23, z=0, nMg=0, dG0=0)
    reacthermo = ReactionThermodynamics(alberty, 'C1')
    reacthermo.pH = 7
    reacthermo.I = 0.1
    reacthermo.T = 298.15
    reacthermo.pMg = 14
    
    c1_reactions = []
    for row in csv.DictReader(open(reaction_fname, 'r')):
        r = Reaction.FromFormula(row['formula'])
        r.Balance(balance_water=False)
        r.SetNames(row['enzyme'])
        dG0_r_prime = float(row['dG0_r_prime'])
        pH, I, T, pMg = [float(row[k]) for k in ['pH', 'I', 'T', 'pMg']]
        reacthermo.AddReaction(r, dG0_r_prime, pH=pH, I=I, T=T, pMg=pMg)
        c1_reactions.append(r)
        
        row['formula'] = r.to_hypertext(show_cids=False)
        dict_list.append(row)
        
    html_writer.write_table(dict_list, headers=['acronym', 'enzyme', 'formula', 'dG0_r_prime'])
    
    reacthermo._Recalculate()
    return reacthermo
    
if __name__ == "__main__":
    html_writer = HtmlWriter('../res/c1_thermodynamics.html')
    reacthermo = GetC1Thermodynamics(html_writer)
    html_writer.close()
    
    formulas = ["C02051 + C00004 => C02972 + C00003",
                "C00143 + C00037 + C00001 => C00101 + C00065",
                "C00288 + 2 C00138 + C00862 => 2 C00001 + 2 C00139 + C01001",
                "C00101 + C00003 => C00415 + C00004"]
    
    reactions = []
    for f in formulas:
        r = Reaction.FromFormula(f)
        r.Balance()
        reactions.append(r)

    print reacthermo.GetTransfromedKeggReactionEnergies(reactions)
