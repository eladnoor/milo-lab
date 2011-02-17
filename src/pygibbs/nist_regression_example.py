from pygibbs.thermodynamics import CsvFileThermodynamics
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist_regression import NistRegression
from toolbox.ods import ODSDictReader
from pygibbs.nist import NistRowData
import pylab

def SingleReaction(sparse):
    html_writer = HtmlWriter("../res/nist/example_reaction.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    nist_regression = NistRegression(db, html_writer, kegg)
    nist_rows = nist_regression.nist.FindRowsAccordingToReaction(sparse)
    dG0_r_tag, ddG0_r, conditions, _, _ = nist_regression.ReverseTranformNistRows(nist_rows)
    dG0_r = dG0_r_tag + ddG0_r
    i = pylab.find(conditions[:, 2] < 9)
    pylab.plot(conditions[i, 2], dG0_r[i, 0], '.')
    pylab.show()

def Pyrophosphatase():
    html_writer = HtmlWriter("../res/nist/example_pyrophosphatase.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    
    nist_rows = []
    for row_dict in ODSDictReader(open('../data/thermodynamics/pyrophosphatase.ods', 'r')):
        row_dict['origin'] = 'pyrophosphatase.ods, line %d' % len(nist_rows)
        if row_dict['comment'].find('Cosolvent = none') == -1:
            continue
        row_dict['comment'] = None
        row_dict['Ktype'] = "Kc'"
        nist_row_data = NistRowData()
        nist_row_data.ReadFromDict(row_dict)
        nist_rows.append(nist_row_data)
    
    html_writer.write("<h2>NIST regression:</h2>")
    nr = NistRegression(db, html_writer, kegg)

    dG0_r_tag, ddG0_r, conditions, _, _ = nr.ReverseTranformNistRows(nist_rows)
    dG0_r = dG0_r_tag + ddG0_r
    pylab.plot(conditions[:, 2], dG0_r, '.')
    pylab.show()
    
if __name__ == "__main__":
    #Pyrophosphatase()
    #SingleReaction({2:-1, 300:-1, 8:1, 2305:1})
    SingleReaction({1:-1, 13:-1, 9:2})