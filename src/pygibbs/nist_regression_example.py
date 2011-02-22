from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist_regression import NistRegression
import pylab
import logging

if __name__ == "__main__":
    html_writer = HtmlWriter("../res/nist/example_reaction.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg.getInstance()
    
    pylab.rcParams['text.usetex'] = False
    pylab.rcParams['legend.fontsize'] = 8
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 8
    pylab.rcParams['lines.linewidth'] = 2
    pylab.rcParams['lines.markersize'] = 5
    pylab.rcParams['figure.figsize'] = [8.0, 6.0]
    pylab.rcParams['figure.dpi'] = 100
    
    nist_regression = NistRegression(db, html_writer)
    
    reactions = {}
    reactions['creatine kinase'] = {2:-1, 300:-1, 8:1, 2305:1}
    reactions['pyrophosphatase'] = {1:-1, 13:-1, 9:2}
    reactions['fructose-bisphosphatase'] = {1:-1, 354:-1, 9:1, 85:1}
    reactions['glutamate => 2-oxoglutarate'] = {1:-1, 3:-1, 25:-1, 4:1, 14:1, 26:1, 80:1}
    reactions['choline-phosphatase'] = {1:-1, 588:-1, 9:1, 114:1}
    reactions['2 ADP => ATP + AMP'] = {8:-2, 2:1, 20:1}
    
    for name, sparse in reactions.iteritems():
        logging.info("Analyzing reaction: " + name)
        html_writer.write('<h2>%s</h2>\n' % name)
        nist_regression.AnalyseSingleReaction(sparse) 
