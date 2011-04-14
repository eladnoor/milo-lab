from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist_regression import NistRegression
from pygibbs.kegg_reaction import Reaction

if __name__ == "__main__":
    html_writer = HtmlWriter("../res/nist/example_reaction.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    
    kegg = Kegg.getInstance()
    nist_regression = NistRegression(db, html_writer)
    nist_regression.nist.T_range = (273.15 + 24, 273.15 + 40)
    #nist_regression.nist.override_I = 0.25
    #nist_regression.nist.override_pMg = 10.0
    
    reactions = []
    #reactions.append(Reaction('methyl-THF hydrolase', {1:-1, 445:-1, 234:1}))
    reactions.append(Reaction('D-glucose kinase', {2:-1, 31:-1, 8:1, 92:1}))
    #reactions.append(Reaction('creatine kinase', {2:-1, 300:-1, 8:1, 2305:1}))
    #reactions.append(Reaction('pyrophosphatase', {1:-1, 13:-1, 9:2}))
    #reactions.append(Reaction('fructose-bisphosphatase', {1:-1, 354:-1, 9:1, 85:1}))
    #reactions.append(Reaction('glutamate => 2-oxoglutarate', {1:-1, 3:-1, 25:-1, 4:1, 14:1, 26:1, 80:1}))
    #reactions.append(Reaction('choline-phosphatase', {1:-1, 588:-1, 9:1, 114:1}))
    #reactions.append(Reaction('2 ADP => ATP + AMP', {8:-2, 2:1, 20:1}))
    #reactions.append(Reaction('galactose dehydrogenase', {124:-1, 3:-1, 4:1, 3383:1}))
    
    for r in reactions:
        r.Balance()
        nist_regression.AnalyzeSingleReaction(r) 
