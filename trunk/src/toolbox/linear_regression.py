import logging
import numpy as np

class LinearRegression(object):
    
    @staticmethod
    def MatrixRank(A, eps=1e-10):
        _U, s, _V = np.linalg.svd(A, full_matrices=True)
        r = len(np.where(s > eps)[0]) # the rank of A
        return r
    
    @staticmethod
    def Kernel(A, reduced_row_echlon=False, eps=1e-10):
        n, m = A.shape
        if n > m:
            return LinearRegression.Kernel(A.T).T
        
        _U, s, V = np.linalg.svd(A, full_matrices=True)
        r = len(np.where(s > eps)[0]) # the rank of A
        
        kerA = np.zeros((m-r, m))
        kerA[0:(m-r), :] = V[r:m,:]

        if reduced_row_echlon:
            LinearRegression.ToReducedRowEchelonForm(kerA)

        return kerA
    
    @staticmethod
    def LeastSquares(A, y, reduced_row_echlon=True, eps=1e-10):
        """
            Performs a safe LeastSquares.
            
            Returns (x, K):
                If A is fully ranked, K = [] and x is the ordinary least squares solution for A*x = y.
                
                If A is not fully ranked, LeastSquares returns the best solution for A*x = y, using the pseudoinverse of cov(A).
                K is the kernel of A, its dimension is: m - rank(A).
                If one wants to interpolate the y-value for a new data point, it must be orthogonal to K (i.e. K*x = 0).
        """
        n, m = A.shape
        if len(y.shape) > 1 and y.shape[1] != 1:
            raise Exception('y is not a column vector')
        if y.shape[0] != n:
            raise Exception('The length of y (%d) does not match the number of rows in A (%d)' % (y.shape[0], n))
        y = np.matrix(y.reshape(n, 1))
        
        zero_columns = np.nonzero([np.linalg.norm(A[:,i]) <= eps for i in xrange(m)])[0]
        nonzero_columns = np.nonzero([np.linalg.norm(A[:,i]) > eps for i in xrange(m)])[0]
        A_red = A[:, nonzero_columns]
        m_red = len(nonzero_columns)
        
        U, s, V = np.linalg.svd(A_red, full_matrices=True)

        r = len(np.where(s > eps)[0]) # the rank of A
        if r < m:
            logging.debug('The rank of A (%d) is lower than the number of columns'
                          ' (%d), i.e. there is a deficiency of dimension %d' % (r, m, m - r))

        inv_S = np.zeros((m_red, n))
        for i in xrange(r):
            inv_S[i, i] = 1/s[i]
        
        x = np.dot(V.T, np.dot(inv_S, np.dot(U.T, y)))
        weights = np.zeros((m, 1))
        weights[nonzero_columns, :] = x

        kerA = np.zeros((m-r, m))
        kerA[0:(m_red-r), nonzero_columns] = V[r:m_red,:]
        for i, j in enumerate(zero_columns):
            kerA[m_red-r+i, j] = 1

        if reduced_row_echlon:
            LinearRegression.ToReducedRowEchelonForm(kerA)

        return weights, kerA
    
    @staticmethod
    def ToReducedRowEchelonForm(A, eps=1e-10, reduce_last_column=True):
        """
            Copied from:
            http://rosettacode.org/wiki/Reduced_row_echelon_form
            
            Slightly changed to fit numpy arrays and an option not to
            reduce the last column (used for solving linear systems)
        """
        lead = 0
        rowCount = A.shape[0]
        columnCount = A.shape[1]
        if not reduce_last_column:
            columnCount -= 1
        
        leads = []
        non_leads = []
        
        for r in xrange(rowCount):
            if lead >= columnCount:
                return leads, non_leads
            i = r
            while abs(A[i, lead]) < eps:
                i += 1
                if i == rowCount:
                    i = r
                    non_leads.append(lead)
                    lead += 1
                    if columnCount == lead:
                        return leads, non_leads
            A[[i,r], :] = A[[r,i], :]  # Replace rows i and r
            A[r, :] /= A[r, lead]      # Divide row r by A[r, lead]
            leads.append(lead)
            for i in xrange(rowCount):
                if i != r and A[i, lead]:
                    A[i, :] -= A[i, lead] * A[r, :] # Subtract for r from row i
            lead += 1
        
        return leads, non_leads
 
    @staticmethod
    def LeastSquaresWithFixedPoints(A, y, index2value):
        """
            The idea is to first solve the Least Squares problem and then
            use the kernel to move the solution to a point where the extra
            constraints are satisfied.
            
            Adding anything of the type ker(A) * a to the solution does not
            change the value of A*x. Therefore, we solve the set of equations:
            
            for each i in index2value: ker(A)*a[i] = index2value[i] - x[i]
            
        """
        x, K = LinearRegression.LeastSquares(A, y, reduced_row_echlon=False)
        K_trunc = K.T[index2value.keys(), :]
        K_inv = np.linalg.pinv(K_trunc)
        if LinearRegression.MatrixRank(K_inv) < len(index2value):
            raise Exception("Cannot set variables to desired values, since the "
                            "truncated kernel matrix is not fully ranked and "
                            "thus not invertible.")
        
        delta_x = np.array([v - x[i] for i,v in index2value.iteritems()])
        a = np.dot(K_inv, delta_x)
        x += np.dot(K.T, a)
        return x, K
 
    @staticmethod
    def SolveLinearSystem(A, b):
        x, _residues, _rank, _s = np.linalg.lstsq(A, b, rcond=1e-10)
        x = x.reshape((A.shape[1], 1)) 
        new_b = np.dot(A, x)
        
        M = np.hstack([A, new_b])
        leads, nonleads = LinearRegression.ToReducedRowEchelonForm(M, reduce_last_column=False)
        
        # since we know the linear system is achievable, all the last rows should be filled with only zeros
        # including the last column (observations).
    
        x = np.zeros((A.shape[1], 1))
        for i, col in enumerate(leads):
            x[col, 0] = M[i, -1]
        
        K = np.zeros((len(nonleads), A.shape[1]))
        row = 0
        for col in xrange(A.shape[1]):
            if col in leads:
                row += 1
            else:
                for i in xrange(row):
                    K[col-row, leads[i]] = -M[i, col].T
                K[col-row, col] = 1
                
        return x, K
    
    class LinearProgrammingException(Exception):
        pass
    
    @staticmethod
    def FindKernel(A):
        b = np.zeros((A.shape[0], 1))
        _x, K = LinearRegression.SolveLinearSystem(A, b)
        return K

if __name__ == '__main__':
#    A = np.matrix([[1, 0, 1, 1],[0, 1, 1, 1],[1, 1, 2, 2]])
#    K = LinearRegression.FindKernel(A)
#    print A
#    print K
#    
#    #y = matrix([[1],[1],[1],[1]])
#    A = np.matrix([[0, 0, 0, 0, 0],[1, 0, 0, 1, 2],[2, 0, 0, 2, 4],[3,0,0,3, 6],[4,0,0,4, 8]])
#    w = np.matrix([[1],[1],[1],[2],[1]])
#
#    w_pred, V = LinearRegression.SolveLinearSystem(A, A*w)
#    print w_pred
#    print V
#    print np.dot(A, V.T)
#    
#    x1 = np.matrix([[0,0,0,1,0]]).T
#    x2 = np.matrix([[-1,0,0,-1,-2]]).T
#    
#    print np.linalg.norm(V*x1)<1e-10, np.dot(x1.T, w_pred)[0,0], np.dot(x1.T, w)[0,0]
#    print np.linalg.norm(V*x2)<1e-10, np.dot(x2.T, w_pred)[0,0], np.dot(x2.T, w)[0,0]
    S = np.matrix([[0,0,1,0,-1,1,0,0],[0,-1,0,-1,1,0,0,0],[-2,0,0,1,0,0,0,0],[0,0,0,0,0,0,-1,1]])
    
    b = np.matrix([[10],[20],[30],[15]])
    f = {0:5, 1:-99, 2:1001, 6:0}
    x, K = LinearRegression.LeastSquaresWithFixedPoints(S, b, f)
    print "solution = ", ', '.join(['%.1f' % i for i in x])
    print "%2s:  %7s  %7s" % ('#', 'fixed', 'sol')
    for key, value in f.iteritems():
        print "%2d:  %7.1f  %7.1f" % (key, value, x[key, 0])
        
    #print K
    
