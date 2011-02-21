from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist_regression import NistRegression
import pylab
import logging

def SingleReaction(nist_regression, sparse):
    logging.info("Analyzing this reaction: " + kegg.sparse_reaction_to_string(sparse))
    nist_rows = nist_regression.nist.FindRowsAccordingToReaction(sparse)
    data = nist_regression.ReverseTranformNistRows(nist_rows)
    
    hyper = kegg.sparse_to_hypertext(sparse, show_cids=False)
    nist_regression.html_writer.write('<h2>Reaction: %s</h2>\n' % hyper)
    fig = pylab.figure()
    nist_regression.html_writer.write('Standard deviations:</br>\n<ul>\n')
    for j, y_axis in enumerate(['dG0_r_tag', 'dG0_r']):
        sigma = pylab.std(data[y_axis])
        nist_regression.html_writer.write("  <li>stdev(%s) = %.2g</li>" % (y_axis, sigma))
        for i, x_axis in enumerate(['pH', 'I', 'pMg']):
            pylab.subplot(2,3,i+3*j+1)
            pylab.plot(data[x_axis], data[y_axis], '.')
            if j == 1:
                pylab.xlabel(x_axis)
            if i == 0:
                pylab.ylabel(y_axis)
    nist_regression.html_writer.write('</ul>\n')
    nist_regression.html_writer.embed_matplotlib_figure(fig, width=640, height=480)
    
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
    SingleReaction(nist_regression, {2:-1, 300:-1, 8:1, 2305:1}) # creatine kinase
    SingleReaction(nist_regression, {1:-1, 13:-1, 9:2}) # pyrophosphatase
    SingleReaction(nist_regression, {1:-1, 354:-1, 9:1, 85:1}) # fructose-bisphosphatase
    SingleReaction(nist_regression, {1:-1, 3:-1, 25:-1, 4:1, 14:1, 26:1, 80:1}) # glutamatse => 2-oxoglutarate
    SingleReaction(nist_regression, {1:-1, 588:-1, 9:1, 114:1}) # choline phosphatase
    SingleReaction(nist_regression, {8:-2, 2:1, 20:1}) # 2 ADP => ATP + AMP