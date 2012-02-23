import logging
import numpy as np
import sys

class LinearRegression(object):
    
    @staticmethod
    def MatrixRank(A, eps=1e-10):
        _U, s, _V = np.linalg.svd(A, full_matrices=False)
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
        m, n = A.shape
        l, k = y.shape
        if l == 1 and k == n: # find a solution for x*A = y (x, y are row vectors)
            return LinearRegression._LeastSquares(A, y, reduced_row_echlon, eps)
        elif k == 1 and l == m: # find a solution for A*x = y (x, y are column vectors)
            weights, kerA = LinearRegression._LeastSquares(A.T, y.T, reduced_row_echlon, eps)
            return weights.T, kerA.T
        raise Exception('The shape of y (%d x %d) does not match the number of '
                        'rows or columns in A (%d x %d)' % (l, k, m, n))
    
    @staticmethod
    def _LeastSquares(A, y, reduced_row_echlon=True, eps=1e-10):
        """
            Performs a safe LeastSquares, by minimizing   || xA - y ||
            
            Returns (x, K):
                If A is fully ranked, K = [] and x is the ordinary least squares solution for xA = y.
                
                If A is not fully ranked, LeastSquares returns the best solution for xA = y, using the pseudoinverse of cov(A).
                K is the kernel of A, its dimension is: n - rank(A).
                If one wants to interpolate the y-value for a new data point, it must be orthogonal to K (i.e. xK = 0).
        """
        U, s, V = np.linalg.svd(A, full_matrices=True)
        r = (s > eps).nonzero()[0].shape[0] # the rank of A
        inv_S = np.matrix(np.diag([1.0/s[i] for i in xrange(r)]))
        V_r = V[:r, :]
        U_r = U[:, :r]
        x = y * V_r.T * inv_S * U_r.T
        kerA = U[:,r:]
        if reduced_row_echlon:
            LinearRegression.ToReducedRowEchelonForm(kerA)
        return x, kerA
    
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
        
        delta_x = np.matrix([v - x[i, 0] for i,v in index2value.iteritems()]).T
        x += K.T * K_inv * delta_x
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
    
    @staticmethod
    def ColumnProjection(A, eps=1e-10):
        """
            Given a matrix A, calculates two projection matrices.
            
            Arguments:
                A - a 2D NumPy array of size n x m
                eps - the eigenvalue threshold for calculating rank (default: 1e-10)
            
            Returns: (P_C, P_L)
                P_C is the projection onto the column-space of A
                P_L is the projection onto the left-null-space of the columns of A
            
            Note:
                P_C + P_L = I
                Both are of size n x n
                rank(P_C) = rank(A)
                rank(P_L) = n - rank(A) 
        """
        P_C, P_L = LinearRegression.RowProjection(A.T, eps)
        return P_C, P_L

    @staticmethod
    def RowProjection(A, eps=1e-10):
        """
            Given a matrix A, calculates two projection matrices.
            
            Arguments:
                A - a 2D NumPy array of size n x m
                eps - the eigenvalue threshold for calculating rank (default: 1e-10)
            
            Returns: (P_R, P_N)
                P_R is the projection onto the row-space of A
                P_N is the projection onto the null-space of the rows of A
            
            Note:
                P_R + P_N = I
                Both are of size m x m
                rank(P_R) = rank(A)
                rank(P_N) = m - rank(A) 
        """
        _U, s, V = np.linalg.svd(A, full_matrices=True)
        r = len(np.where(s > eps)[0]) # the rank of A
        P_R = V[:r,:].T * V[:r,:] # a projection matrix onto the row-space of A
        P_N = V[r:,:].T * V[r:,:] # a projection matrix onto the null-space of A
        return P_R, P_N

    @staticmethod
    def ColumnUniqueRegression(A, y):
        n = A.shape[1]
        assert np.prod(y.shape) == n
        y = y.reshape((1, n))

        A_unique, col_mapping = LinearRegression.ColumnUnique(A)
        y_unique = np.zeros((1, len(col_mapping)))

        for i, col_indices in enumerate(col_mapping):
            y_unique[0, i] = np.mean(y[0, col_indices])
        
        return A_unique, y_unique

    @staticmethod
    def RowUniqueRegression(A, y):
        """
            A procedure usually performed before linear regression (i.e. solving Ax = y).
            If the matrix A contains repeating rows, it is advisable to combine
            all of them to one row, and the observed value corresponding to that
            row will be the average of the original observations.
            
            Returns:
                A_unique, y_unique
                
                where A_unique has the same number of columns as A, but with
                unique rows, and y_unique is shorter as well.
        """
        n = A.shape[0]
        assert np.prod(y.shape) == n
        y = y.reshape((n, 1))

        A_unique, row_mapping = LinearRegression.RowUnique(A)
        y_unique = np.zeros((len(row_mapping), 1))

        for i, row_indices in enumerate(row_mapping):
            y_unique[i, 0] = np.mean(y[row_indices, 0])
        
        return A_unique, y_unique

    @staticmethod
    def RowUnique(A, remove_zero=False):
        """
            A procedure usually performed before linear regression (i.e. solving Ax = y).
            If the matrix A contains repeating rows, it is advisable to combine
            all of them to one row, and the observed value corresponding to that
            row will be the average of the original observations.

            Input:
                A - a 2D NumPy array
                remove_zero - if True, removes any only-zero rows from the result
            
            Returns:
                A_unique, row_mapping
                
                where A_unique has the same number of columns as A, but with
                unique rows.
                row_mapping is a dictionary whose keys are rows in A_unique and
                whose values are lists of corresponding rows in A.
        """
        A_unique = np.unique([tuple(A[i,:].flat)
                              for i in xrange(A.shape[0])
                              if (not remove_zero or (A[i,:] != 0).any())])
        A_unique = np.array(A_unique)

        row_mapping = {}
        for i in xrange(A_unique.shape[0]):
            # find the indices of rows in A which correspond to this unique row in A_unique
            row_vector = A_unique[i:i+1,:]
            diff = abs(A - np.repeat(row_vector, A.shape[0], 0))
            row_mapping[i] = list(np.where(np.sum(diff, 1) == 0)[0].flat)
        
        return A_unique, row_mapping
    
    @staticmethod
    def ColumnUnique(A, remove_zero=False):
        A_unique, column_mapping = LinearRegression.RowUnique(A.T, remove_zero)
        return A_unique.T, column_mapping
    
if __name__ == '__main__':
    eps = 1e-10
    A = np.matrix([[1,0,0],[0,1,0],[0,0,0],[1,0,0],[0,1,0],[0,1,0]])
    print LinearRegression.RowUnique(A, remove_zero=True)
    b = np.matrix([[1],[2],[1],[3],[4],[1]])
    #A, b = LinearRegression.RowUniqueRegression(A, b)
    print 'A:\n', A
    print 'b:\n', b.T.tolist()[0]
    
    G = np.matrix([[1],[2],[0]])
    P_C, P_L = LinearRegression.ColumnProjection(A)
    print "P_C:\n", P_C
    print "P_C * P_C^T:\n", np.dot(P_C, P_C.T)
    
    P_C2, P_L2 = LinearRegression.ColumnProjection(A*G)

    n = A.shape[0]
    print "b:\n", b.T.tolist()[0]
    b_C = np.dot(P_C, b)
    b_L = np.dot(P_L, b)
    print "b_C:\n", b_C.T.tolist()[0]
    print "b_L:\n", b_L.T.tolist()[0]

    b_C2 = np.dot(P_C2, b)
    b_L2 = np.dot(P_L2, b)
    print "b_C2:\n", b_C2.T.tolist()[0]
    print "b_L2:\n", b_L2.T.tolist()[0]

    #x, _ = LinearRegression.LeastSquares(A, b)
    #print "A*x - b:", (np.dot(A, x) - b).T.tolist()[0]

    x_C, _ = LinearRegression.LeastSquares(A*G, b_C)
    print "A*x_C - b_C:\n", abs(A * G * x_C - b_C).T.tolist()[0]
    print abs((P_C - P_C2) * b).T.tolist()[0]
    
    #x_N, _ = LinearRegression.LeastSquares(A, b_N)
    #print "A*x_N - b_N:", (np.dot(A, x_N) - b_N).T.tolist()[0]

    #print (P_R * b).T.tolist()[0]
    #print (P_R2 * b).T.tolist()[0]
    #print (P_R2 * P_R * b).T.tolist()[0]
    #print (P_R * P_R2 * b).T.tolist()[0]
    #sys.exit(0)

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
    f = {0:5, 1:-99, 2:1021, 6:0}
    x, K = LinearRegression.LeastSquaresWithFixedPoints(S, b, f)
    print "solution = ", ', '.join(['%.1f' % i for i in x])
    print "%2s:  %7s  %7s" % ('#', 'fixed', 'sol')
    for key, value in f.iteritems():
        print "%2d:  %7.1f  %7.1f" % (key, value, x[key, 0])
        
    #print K
    
