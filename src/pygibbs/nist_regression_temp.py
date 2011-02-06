from toolbox.linear_regression import LinearRegression
import pylab
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase

def main():
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    cids = pylab.np.loadtxt('../res/nist/regress_CID.txt', delimiter=',')
    S = pylab.np.loadtxt('../res/nist/regress_S.txt', delimiter=',')
    dG0 = pylab.np.loadtxt('../res/nist/regress_dG0.txt', delimiter=',')

    print S.shape
    
    nonzero_columns = pylab.find(pylab.sum(abs(S), 0))
    S = S[:, nonzero_columns]
    cids = [cids[i] for i in nonzero_columns]
        
    dG0 = dG0.reshape((dG0.shape[0], 1))
    x1, K1 = LinearRegression.LeastSquares(S, dG0)
    x2, K2 = LinearRegression.SolveLinearSystem(S, dG0)

    print pylab.norm(pylab.dot(S, K1.T))
    print pylab.norm(pylab.dot(S, K2.T))
        
    for i in xrange(K2.shape[0]):
        nonzero_columns = pylab.find(abs(K2[i, :]) > 1e-10)
        gv = " + ".join(["%g %s (C%05d)" % (K2[i, j], kegg.cid2name(int(cids[j])), cids[j]) for j in nonzero_columns])
        print i, ":", gv
                
if __name__ == "__main__":
    main()