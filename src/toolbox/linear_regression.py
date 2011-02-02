import logging
from numpy.linalg import svd, norm
from numpy import matrix, zeros, dot
from pylab import find

class LinearRegression(object):
    
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
        y = matrix(y.reshape(n, 1))
        
        zero_columns = find([norm(A[:,i])<=eps for i in xrange(m)])
        nonzero_columns = find([norm(A[:,i])>eps for i in xrange(m)])
        A_red = A[:, nonzero_columns]
        m_red = len(nonzero_columns)
        
        U, s, V = svd(A_red, full_matrices=True)

        r = len(find(s > eps)) # the rank of A
        if r < m:
            logging.debug('The rank of A (%d) is lower than the number of columns'
                          ' (%d), i.e. there is a deficiency of dimension %d' % (r, m, m - r))

        inv_S = zeros((m_red, n))
        for i in xrange(r):
            inv_S[i, i] = 1/s[i]
        
        x = dot(V.T, dot(inv_S, dot(U.T, y)))
        weights = zeros((m, 1))
        weights[nonzero_columns, :] = x

        kerA = zeros((m-r, m))
        kerA[0:(m_red-r), nonzero_columns] = V[r:m_red,:]
        for i, j in enumerate(zero_columns):
            kerA[m_red-r+i, j] = 1

        if reduced_row_echlon:
            LinearRegression.GaussJordan(kerA, eps)
            #LinearRegression.ToReducedRowEchelonForm(kerA)

        return weights, kerA
    
    @staticmethod
    def ToReducedRowEchelonForm(M, reduce_last_column=True):
        """
            Copied from:
            http://rosettacode.org/wiki/Reduced_row_echelon_form
            
            Slightly changed to fit numpy arrays and an option not to
            reduce the last column (used for solving linear systems)
        """
        lead = 0
        rowCount = M.shape[0]
        columnCount = M.shape[1]
        if not reduce_last_column:
            columnCount -= 1
        
        for r in xrange(rowCount):
            if lead >= columnCount:
                return
            i = r
            while M[i, lead] == 0:
                i += 1
                if i == rowCount:
                    i = r
                    lead += 1
                    if columnCount == lead:
                        return
            M[[i,r], :] = A[[r,i], :]  # Replace rows i and r
            M[r, :] /= M[r, lead]      # Divide row r by M[r, lead]
            for i in xrange(rowCount):
                if i != r and M[i, lead]:
                    M[i, :] -= M[i, lead] * M[r, :] # Subtract for r from row i
            lead += 1
    
    @staticmethod
    def GaussJordan(A, eps=1e-10):
        """
            Puts given matrix (2D array) into the Reduced Row Echelon Form.
            Returns True if successful, False if 'm' is singular.
            NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
            Written by Jarno Elonen in April 2005, released into Public Domain
        """
        h, w = A.shape
        for y in xrange(0, h):
            maxrow = y
            for y2 in xrange(y+1, h):    # Find max pivot
                if abs(A[y2, y]) > abs(A[maxrow, y]):
                    maxrow = y2
            A[[y,maxrow], :] = A[[maxrow,y], :]
            if abs(A[y, y]) <= eps:     # Singular?
                return False
            for y2 in xrange(y+1, h):    # Eliminate column y
                c = A[y2, y] / A[y, y]
                for x in xrange(y, w):
                    A[y2, x] -= A[y, x] * c
        for y in xrange(h-1, 0-1, -1): # Backsubstitute
            c  = A[y, y]
            for y2 in xrange(0, y):
                for x in xrange(w-1, y-1, -1):
                    A[y2, x] -=  A[y, x] * A[y2, y] / c
            A[y, y] /= c
            for x in xrange(h, w):       # Normalize row y
                A[y, x] /= c
        return True
    
    @staticmethod
    def Rank(A, eps=1e-10):
        _U, s, _V = svd(A, full_matrices=False)
        return len(find(s > eps))
    
if __name__ == '__main__':
    #A = matrix([[1, 2, 3, 9],[2, -1, 1, 8],[2, 4, 6, 3]])
    #LinearRegression.ToReducedRowEchelonForm(A, reduce_last_column=False)
    #print A
    #y = matrix([[1],[1],[1],[1]])
    A = matrix([[0, 0, 0, 0, 0],[1, 0, 0, 1, 2],[2, 0, 0, 2, 4],[3,0,0,3, 6],[4,0,0,4, 8]])
    
    w = matrix([[1],[1],[1],[2],[1]])
    
    #y = A*x
    #print A
    #print x
    w_pred, V = LinearRegression.LeastSquares(A, A*w)
    print w_pred
    print V
    print dot(A, V.T)

    
    x1 = matrix([[0,0,0,1,0]]).T
    x2 = matrix([[-1,0,0,-1,-2]]).T
    
    print norm(V*x1)<1e-10, dot(x1.T, w_pred)[0,0], dot(x1.T, w)[0,0]
    print norm(V*x2)<1e-10, dot(x2.T, w_pred)[0,0], dot(x2.T, w)[0,0]