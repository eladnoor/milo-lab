import logging
from numpy.linalg import svd, norm, inv
from numpy import matrix, zeros, dot, float64, ones
from pylab import find
from sympy import Matrix

class LinearRegression(object):
    
    @staticmethod
    def LeastSquares(A, y, prior=None, reduced_row_echlon=True, eps=1e-10):
        """
            Input:
                A - a regression matrix, of dimension n x m
                y - the solution vector, of dimension n x 1 

            Method:
                Performs a safe LeastSquares, by using SVD to calculate the pseudoinverse of A:
                A = U*D*V   =>   x = V' * D^-1 * U' * y
                where U'*U = I and V'*V = I and D is diagonal (although doesn't have to be square).
            
            Returns (x, K):
                x - the ordinary least squares solution, of dimension m x 1
                
                If rank(A) = m the return value is the unique LeastSquares solution
                If rank(A) < m the return value is one of the LeastSquares solution, where
                    the ker(A) dimensions are arbitrarily assigned a weight of 0.
                
                K = ker(A), where dim(K) = m - rank(A).
            
            Note:
                If one wants to interpolate the y-value for a new data point,
                    it must be orthogonal to K (i.e. K*x = 0).
        """
        n, m = A.shape
        if len(y.shape) > 1 and y.shape[1] != 1:
            raise Exception('y is not a column vector')
        if y.shape[0] != n:
            raise Exception('The length of y (%d) does not match the number of rows in A (%d)' % (y.shape[0], n))
        if prior != None:
            if prior.shape[0] != m and prior[1] != 1:
                raise Exception('The shape of the prior does not match A: %s' % str(prior.shape))
        else:
            prior = zeros((m, 1))
        
        U, s, V = svd(A, full_matrices=True)

        r = len(find(s > eps)) # the rank of A
        if r < m:
            logging.debug('The rank of A (%d) is lower than the number of columns'
                          ' (%d), i.e. there is a deficiency of dimension %d' % (r, m, m - r))

        inv_S = zeros((A.shape[1], A.shape[0]))
        for i in xrange(r):
            inv_S[i,i] = 1/s[i]
        
        x = dot(V.T, dot(inv_S, dot(U.T, y)))
        K = V[r:,:] # an orthonormal basis for ker(A)

        # A normalization that will ensure that zero-columns will have 0 values
        # in the solution
        x = x + dot(dot(K.T, K), -x)

        if reduced_row_echlon:
            #LinearRegression.GaussJordan(K, eps)
            K = LinearRegression.ReducedRowEchelon(K)

        return x, K
    
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
    A = matrix([[1, 0, -1],[2, 0, -2]])
    y = matrix([[1],[1]])
    w_pred, V = LinearRegression.LeastSquares(A, y)
    print w_pred
    print V
    print "-"*100
    
    A = matrix([[0, 0, 0, 0, 0],[1, 0, 0, 1, 2],[2, 0, 0, 2, 4],[3,0,0,3, 6],[4,0,0,4, 8]])
    w = matrix([[1],[1],[1],[2],[1]])
    
    #y = A*x
    #print A
    #print x
    w_pred, V = LinearRegression.LeastSquares(A, A*w, reduced_row_echlon=False)
    print w_pred
    print V
    
    x1 = matrix([[0,0,0,1,0]]).T
    x2 = matrix([[-1,0,0,-1,-2]]).T
    
    print norm(V*x1)<1e-10, dot(x1.T, w_pred)[0,0], dot(x1.T, w)[0,0]
    print norm(V*x2)<1e-10, dot(x2.T, w_pred)[0,0], dot(x2.T, w)[0,0]