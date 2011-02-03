from toolbox.linear_regression import LinearRegression
import pylab

def SolveLinearSystem(A, b):
    x, _residues, _rank, _s = pylab.np.linalg.lstsq(A, b, rcond=1e-10)
    new_b = pylab.dot(A, x)
    
    M = pylab.hstack([A, new_b])
    leads, nonleads = LinearRegression.ToReducedRowEchelonForm(M, reduce_last_column=False)
    
    # since we know the linear system is achievable, all the last rows should be filled with only zeros
    # including the last column (observations).

    x = pylab.zeros((A.shape[1], 1))
    for i, col in enumerate(leads):
        x[col, 0] = M[i, -1]
    
    K = pylab.zeros((len(nonleads), A.shape[1]))
    row = 0
    for col in xrange(A.shape[1]):
        if col in leads:
            row += 1
        else:
            K[col-row, 0:row] = M[0:row, col].T
            K[col-row, col] = 1

    return x, K


def main():
    cids = pylab.np.loadtxt('../res/nist/regress_CID.txt', delimiter=',')
    S = pylab.np.loadtxt('../res/nist/regress_S.txt', delimiter=',')
    dG0 = pylab.np.loadtxt('../res/nist/regress_dG0.txt', delimiter=',')

    print S.shape
        
    dG0 = dG0.reshape((dG0.shape[0], 1))
    x1, K1 = LinearRegression.LeastSquares(S, dG0)
    x2, K2 = SolveLinearSystem(S, dG0)

    for row in xrange(K2.shape[0]):
        print row,
        for col in xrange(K2.shape[1]):
            if abs(K2[row, col]) > 1e-10:
                print '(%d, %d) %g, ' % (row, col, K2[row, col]),
        print ''
                
    #pylab.figure()
    #pylab.plot(pylab.dot(S, x1), dG0, '.')
    #pylab.figure()
    #pylab.plot(pylab.dot(S, x2), dG0, '.')
    #pylab.show()
    

if __name__ == "__main__":
    main()