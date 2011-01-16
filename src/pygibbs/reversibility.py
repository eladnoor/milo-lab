from toolbox import database
from toolbox.html_writer import HtmlWriter
from pygibbs import thermodynamics
from pygibbs.thermodynamic_constants import R, default_I, default_pH,\
    default_pMg, default_T
import pylab
from pygibbs.groups import GroupContribution
import sys
from pygibbs.kegg import Kegg

def CalculateReversability(rid, G, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T):
    dG0 = G.estimate_dG_keggrid(rid, pH, pMg, I, T)
    sparse = G.kegg().rid2sparse_reaction(rid)
    sum_s = sum(sparse.values())
    sum_abs_s = sum([abs(x) for x in sparse.values()])
    
    return 2 * ((-dG0/(R*T) + pylab.log(c_mid)*sum_s) / sum_abs_s)

def main():
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/dG0_test.html')
    kegg = Kegg(db)
    G = GroupContribution(db, html_writer=html_writer, kegg=kegg)
    G.init()
    #kegg = G.kegg()
    #kegg.db = db
    #kegg.ToDatabase()
    #sys.exit(0)
    
    c_mid = 1e-3
    pH, pMg, I, T = (7, 3, 0.1, 298.15)
    
    histogram = {}
    total_hist = []
    hits = 0
    misses = 0
    for rid_flux_list in kegg.mid2rid_map.itervalues():
        if not rid_flux_list or len(rid_flux_list) < 2:
            continue
        for i, (rid, flux) in enumerate(rid_flux_list):
            try:
                r = flux * CalculateReversability(rid, G, c_mid, pH, pMg, I, T)
                histogram.setdefault(i, []).append(r)
                if i > 1:
                    total_hist.append(r)
                hits += 1
            except thermodynamics.MissingCompoundFormationEnergy:
                misses += 1
                continue
    
    print "Reactions with known dG0", hits
    print "Reactions with unknown dG0", misses
    
    medians = []
    for i in histogram.keys():
        if i < 8:
            medians.append(pylab.median(histogram[i]))  
    
    pylab.figure()
    pylab.bar(range(len(medians)), medians)
    pylab.xlabel('Position in pathway')
    pylab.ylabel('Median reversibility')

    pylab.figure()
    pylab.subplot(3,1,1)
    pylab.hist(histogram[0], bins=50, range=(-60, 60))
    pylab.subplot(3,1,2)
    pylab.hist(histogram[1], bins=50, range=(-60, 60))
    pylab.subplot(3,1,3)
    pylab.hist(total_hist, bins=50, range=(-60, 60))
    pylab.show()

if __name__ == "__main__":
    main()