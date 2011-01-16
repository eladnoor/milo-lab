import logging
from numpy.linalg import svd, inv, norm
from numpy import diag, matrix, hstack, zeros, dot, vstack, float64
from pylab import find
from sympy import Matrix
import sys

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
        if n < m:
            raise Exception('Linear problem is under-determined (more variables than equations)')
        
        zero_columns = find([norm(A[:,i])<=eps for i in xrange(m)])
        nonzero_columns = find([norm(A[:,i])>eps for i in xrange(m)])
        A_red = A[:, nonzero_columns]
        m_red = len(nonzero_columns)
        
        U, s, V = svd(A_red, full_matrices=True)

        r = len(find(s > eps)) # the rank of A
        if r < m:
            logging.debug('The rank of A (%d) is lower than the number of columns'
                          ' (%d), i.e. there is a deficiency of dimension %d' % (r, m, m - r))

        inv_V = inv(V)
        inv_U = inv(U)
        
        D = diag([1/s[i] for i in xrange(r)] + [0] * (m_red-r))
        inv_S = hstack([D, zeros((m_red, n-m_red))])
        
        x = dot(inv_V, dot(inv_S, dot(inv_U, y)))
        weights = zeros((m, 1))
        weights[nonzero_columns, :] = x

        kerA = zeros((m-r, m))
        kerA[0:(m_red-r), nonzero_columns] = V[r:m_red,:]
        for i, j in enumerate(zero_columns):
            kerA[m_red-r+i, j] = 1

        if reduced_row_echlon:
            #LinearRegression.GaussJordan(kerA, eps)
            kerA = LinearRegression.ReducedRowEchelon(kerA)

        return weights, kerA
    
    @staticmethod
    def ReducedRowEchelon(A):
        R = Matrix(A).rref()
        return float64(matrix(R[0]))
    
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
    
if __name__ == '__main__':
    #A = matrix([[1, 2, 3],[2, 3, 4],[-1, 8, 2],[2, 3, 1]])
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