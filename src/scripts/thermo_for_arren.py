from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase
from pygibbs.kegg import Kegg
import pylab
from toolbox import util
import csv


def main():
    pH, pMg, I, T = (7.0, 3, 0.1, 298.15)
    
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg.getInstance()
    alberty = PsuedoisomerTableThermodynamics('../data/thermodynamics/alberty_pseudoisomers.csv')
    
    cids = alberty.get_all_cids()
    dG0_f = pylab.zeros((len(cids), 1))

    for i, cid in enumerate(cids):
        dG0_f[i, 0] = alberty.cid2dG0_tag(cid, pH=pH, pMg=pMg, I=I, T=T)
    
    S = pylab.zeros((0, len(cids)))
    rids = []
    ec_numbers = []
    
    for rid in kegg.get_all_rids():
        sparse = kegg.rid2sparse_reaction(rid)
        if not set(cids).issuperset(sparse.keys()):
            continue
        
        rids.append(rid)
        ec_numbers.append(kegg.rid2ec_list(rid))
        S_row = pylab.zeros((1, len(cids)))
        for cid, coeff in sparse.iteritems():
            S_row[0, cids.index(cid)] = coeff
        S = pylab.vstack([S, S_row])
    
    dG0_r = pylab.dot(S, dG0_f)

    util._mkdir('../res/arren')
    s_writer = csv.writer(open('../res/arren/stoichiomety.csv', 'w'))
    r_writer = csv.writer(open('../res/arren/reactions.csv', 'w'))
    e_writer = csv.writer(open('../res/arren/ec_numbers.csv', 'w'))
    r_writer.writerow(['rid', 'dG0_r'])
    e_writer.writerow(['rid', 'ec0', 'ec1', 'ec2', 'ec3'])
    for i in xrange(S.shape[0]):
        s_writer.writerow(["%d" % x for x in S[i,:]])
        for ec in ec_numbers[i].split(';'):
            e_writer.writerow(['%d' % rids[i]] + ec.split('.'))
        r_writer.writerow(["%d" % rids[i], '%.1f' % dG0_r[i,0]])
    
    c_writer = csv.writer(open('../res/arren/compounds.csv', 'w'))
    c_writer.writerow(['cid', 'dG0_f'])
    for j in xrange(len(cids)):
        c_writer.writerow(['%d' % cids[j], '%.1f' % dG0_f[j, 0]])
    
if __name__ == "__main__":
    main()