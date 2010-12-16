import logging
from numpy.linalg import svd, inv
from numpy import diag, matrix, hstack, zeros, array, dot
from pylab import find

class LinearRegression:
    
    @staticmethod
    def LeastSquares(A, y):
        """
            Performs a safe LeastSquars.
            
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
        
        U, s, V = svd(A, full_matrices=True)

        r = len(find(s > 1e-8)) # the rank of A
        if r < m:
            logging.warning('The rank of A (%d) is lower than the number of columns'
                            ' (%d), i.e. there is a deficiency of dimension %d' % (r, m, m - r))

        inv_V = inv(V)
        inv_U = inv(U)
        
        D = diag([1/s[i] for i in xrange(r)] + [0] * (m-r))
        inv_S = hstack([D, zeros((m, n-m))])
        
        print inv_V.shape, inv_S.shape, inv_U.shape, y.shape
        
        x = dot(inv_V, dot(inv_S, dot(inv_U, y)))
        return x, V[r:m,:]
    
    
if __name__ == '__main__':
    #A = matrix([[1, 2, 3],[2, 3, 4],[-1, 8, 2],[2, 3, 1]])
    #y = matrix([[1],[1],[1],[1]])
    A = matrix([[0, 1, 1],[1, 2, 2],[-1, 1, 1]])
    w = matrix([[1],[1],[1]])
    y = array([2, 5, 1])
    
    #y = A*x
    #print A
    #print x
    w_pred, V = LinearRegression.LeastSquares(A, A*w)
    print w_pred
    
    x1 = matrix([[0],[5],[5]])
    x2 = matrix([[0],[1],[0]])
    
    print V*x1, x1.T*w_pred
    print V*x2, x2.T*w_pred