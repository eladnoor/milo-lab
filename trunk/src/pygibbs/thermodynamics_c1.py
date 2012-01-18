from pygibbs.thermodynamics import ReactionThermodynamics,\
    PsuedoisomerTableThermodynamics, BinaryThermodynamics
import csv
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg import Kegg
from toolbox.html_writer import HtmlWriter
from toolbox.database import SqliteDatabase

def GetC1Thermodynamics(html_writer, reaction_fname='../data/thermodynamics/c1_reaction_thermodynamics.csv'):
    #kegg = Kegg.getInstance()
    
    html_writer.write("<h1>C1 thermodynamics</h1>\n")
    
    dict_list = []
    db_public = SqliteDatabase('../data/public_data.sqlite')
    alberty = PsuedoisomerTableThermodynamics.FromDatabase(\
                        db_public, 'alberty_pseudoisomers', name='alberty')
    reacthermo = ReactionThermodynamics(alberty, 'C1')
    reacthermo.pH = 7
    reacthermo.I = 0.1
    reacthermo.T = 298.15

    c1_reactions = []
    for row in csv.DictReader(open(reaction_fname, 'r')):
        r = Reaction.FromFormula(row['formula'])
        r.Balance(balance_water=False)
        r.SetNames(row['enzyme'])
        dG0_r_prime = float(row['dG0_r_prime'])
        reacthermo.AddReaction(r, dG0_r_prime)
        c1_reactions.append(r)
        
        row['formula'] = r.to_hypertext(show_cids=False)
        dict_list.append(row)
        
    html_writer.write_table(dict_list, headers=['acronym', 'enzyme', 'formula', 'dG0_r_prime'])
    
    reacthermo._Recalculate()
    print reacthermo.GetTransfromedKeggReactionEnergies(c1_reactions)
    
    return reacthermo
    
if __name__ == "__main__":
    html_writer = HtmlWriter('../res/c1_thermodynamics.html')
    reacthermo = GetC1Thermodynamics(html_writer)
    html_writer.close()
    
    r = Reaction.FromFormula("2 C00143 => 2 C00445")
    r.Balance()
    print reacthermo.GetTransfromedKeggReactionEnergies([r])
