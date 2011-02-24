from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.nist_regression import NistRegression
from pygibbs.dissociation_constants import DissociationConstants
import logging

if __name__ == "__main__":
    html_writer = HtmlWriter("../res/nist/example_reaction.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    
    dissociation = DissociationConstants(db, html_writer)
    dissociation.LoadValuesToDB()

    kegg = Kegg.getInstance()
    nist_regression = NistRegression(db, html_writer)
    T_range = (298, 314)
    
    reactions = {}
    #reactions['creatine kinase'] = {2:-1, 300:-1, 8:1, 2305:1}
    reactions['pyrophosphatase'] = {1:-1, 13:-1, 9:2}
    #reactions['fructose-bisphosphatase'] = {1:-1, 354:-1, 9:1, 85:1}
    #reactions['glutamate => 2-oxoglutarate'] = {1:-1, 3:-1, 25:-1, 4:1, 14:1, 26:1, 80:1}
    #reactions['choline-phosphatase'] = {1:-1, 588:-1, 9:1, 114:1}
    #reactions['2 ADP => ATP + AMP'] = {8:-2, 2:1, 20:1}
    #reactions['galactose dehydrogenase'] = {124:-1, 3:-1, 4:1, 3383:1}
    
    for name, sparse in reactions.iteritems():
        logging.info("Analyzing reaction: " + name)
        html_writer.write('<h2>%s</h2>\n' % name)
        
        nist_rows = nist_regression.nist.SelectRowsFromNist(sparse, T_range)
        dict_list = []
        for nist_row_data in nist_rows:
            d = {}
            d['pH'] = nist_row_data.pH
            d['I'] = nist_row_data.I
            d['pMg'] = nist_row_data.pMg
            d['dG\'0_r'] = "%.2f" % nist_row_data.dG0_r
            d['T(K)'] = nist_row_data.T
            dict_list.append(d)
        html_writer.write_table(dict_list, headers=['T(K)', 'pH', 'I', 'pMg', 'dG\'0_r'])
        html_writer.write('</br>\n')
        
        nist_regression.AnalyseSingleReaction(sparse, T_range) 
        for cid, coeff in sparse.iteritems():
            for pseudoisomer in nist_regression.cid2diss_table[cid].GenerateAll():
                print "C%05d" % cid, pseudoisomer