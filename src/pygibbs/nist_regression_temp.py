from toolbox.linear_regression import LinearRegression
from pygibbs.kegg import Kegg
import pylab

def main():
    kegg = Kegg.getInstance()
    cids = pylab.np.loadtxt('../res/nist/regress_CID.txt', delimiter=',')
    S = pylab.np.loadtxt('../res/nist/regress_S.txt', delimiter=',')
    dG0 = pylab.np.loadtxt('../res/nist/regress_dG0.txt', delimiter=',')

    print "S has %d rows and %d columns" % S.shape
    
    nonzero_columns = pylab.find(pylab.sum(abs(S), 0))
    S = S[:, nonzero_columns]
    cids = [cids[i] for i in nonzero_columns]
        
    dG0 = dG0.reshape((dG0.shape[0], 1))
    #x, K = LinearRegression.SolveLinearSystem(S, dG0)
    try:
        K = LinearRegression.FindSparseKernel(S)
    except LinearRegression.LinearProgrammingException as e:
        print "Error when trying to find a sparse kernel: " + str(e)

    for i in xrange(K.shape[0]):
        nonzero_columns = pylab.find(abs(K[i, :]) > 1e-10)
        gv = " + ".join(["%g %s (C%05d)" % (K[i, j], kegg.cid2name(int(cids[j])), cids[j]) for j in nonzero_columns])
        print i, ":", gv
                
if __name__ == "__main__":
    main()