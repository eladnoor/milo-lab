from toolbox.database import SqliteDatabase
from pygibbs.kegg import Kegg
from toolbox.html_writer import HtmlWriter
from pygibbs.groups import GroupContribution
from pygibbs.nist import Nist
import logging

def main():    
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    html_writer = HtmlWriter('../res/groups_quicktest.html')
    
    nist = Nist(db, html_writer, kegg)
    nist.Load()

    G = GroupContribution(db=db, html_writer=html_writer, kegg=kegg)
    G.load_groups("../data/thermodynamics/groups_species.csv")
    if True:
        G.train("../data/thermodynamics/dG0.csv", use_dG0_format=True)
    else:
        G.load_training_data()
    G.override_gc_with_measurements = True
    G.quick_init(nist.GetAllCids(), html_fname='../res/groups_quicktest_pmaps.html')

    logging.info('Writing the NIST report')
    num_estimations, rmse = nist.verify_results(G, T_range=(298, 314))
    logging.info('N = %d' % num_estimations)
    logging.info('RMSE = %.1f' % rmse)

if __name__ == '__main__':
    main()