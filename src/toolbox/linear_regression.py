import logging
from pylab import np, find, array

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
        y = np.matrix(y.reshape(n, 1))
        
        zero_columns = find([np.linalg.norm(A[:,i])<=eps for i in xrange(m)])
        nonzero_columns = find([np.linalg.norm(A[:,i])>eps for i in xrange(m)])
        A_red = A[:, nonzero_columns]
        m_red = len(nonzero_columns)
        
        U, s, V = np.linalg.svd(A_red, full_matrices=True)

        r = len(find(s > eps)) # the rank of A
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
    def Rank(A, eps=1e-10):
        _U, s, _V = np.linalg.svd(A, full_matrices=False)
        return len(find(s > eps))

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
    def FindSparseKernel(A, kernel_dimension=None, eps=1e-10, upper_bound=1000):
        """
            Finds a sparse representation of the kernel matrix, using MILP
            to iterate the Fundamental Modes of the matrix.
        
            Input:
                a (n x m) matrix A, whose rank is r.
            Return:
                a (m x m-r) matrix K that will span the kernel of A, i.e.:
                span(K) = {x | Ax = 0}
        """
        import cplex
        
        kernel_rank = A.shape[1] - LinearRegression.Rank(A, eps)
        if not kernel_dimension or kernel_dimension > kernel_rank:
            kernel_dimension = kernel_rank
        cpl = cplex.Cplex()
        cpl.set_problem_name('find_kernel')
        cpl.set_log_stream(None)
        cpl.set_results_stream(None)
        cpl.set_warning_stream(None)
        
        for col in xrange(A.shape[1]):
            cpl.variables.add(names=['c%d_plus' % col], lb=[0], ub=[upper_bound])
            cpl.variables.add(names=['g%d_plus' % col], types='B')

            # c_plus - M*g_plus <= 0
            constraint_name = 'g%d_plus_bound_upper' % col
            cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_plus' % col, 1)
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, -upper_bound)

            # c_plus - g_plus >= 0
            constraint_name = 'g%d_plus_bound_lower' % col
            cpl.linear_constraints.add(names=[constraint_name], senses='G', rhs=[0])
            cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_plus' % col, 1)
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, -1)

            cpl.variables.add(names=['c%d_minus' % col], lb=[0], ub=[upper_bound])
            cpl.variables.add(names=['g%d_minus' % col], types='B')

            # c_minus - M*g_minus <= 0
            constraint_name = 'g%d_minus_bound_upper' % col
            cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_minus' % col, 1)
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, -upper_bound)

            # c_minus - g_minus >= 0
            constraint_name = 'g%d_minus_bound_lower' % col
            cpl.linear_constraints.add(names=[constraint_name], senses='G', rhs=[0])
            cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_minus' % col, 1)
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, -1)
            
            # g_plus + g_minus <= 1
            constraint_name = 'g%d_bound' % col
            cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[1])
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, 1)
            cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, 1)

        for row in xrange(A.shape[0]):
            constraint_name = 'r%d' % row
            cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
            for col in xrange(A.shape[1]):
                if A[row,col] != 0:
                    cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_plus' % col, A[row,col])
                    cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_minus' % col, -A[row,col])

        all_gammas = ['g%d_plus' % col for col in xrange(A.shape[1])] + \
                     ['g%d_minus' % col for col in xrange(A.shape[1])]
        
        cpl.objective.set_linear([(g, 1) for g in all_gammas])

        cpl.linear_constraints.add(names=['avoid_0'], senses='G', rhs=[1])
        for g in all_gammas:
            cpl.linear_constraints.set_coefficients('avoid_0', g, 1)

        K = np.zeros((A.shape[1], kernel_dimension))
        dimension = 0
        emf_counter = 0
        while dimension < kernel_dimension:
            emf_counter += 1
            try:
                cpl.solve()
            except cplex.exceptions.CplexSolverError as e:
                raise LinearRegression.LinearProgrammingException(str(e))
            
            if cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.MIP_optimal:
                raise LinearRegression.LinearProgrammingException("No more EMFs")
            
            g_plus = cpl.solution.get_values(['g%d_plus' % col for col in xrange(A.shape[1])])
            g_minus = cpl.solution.get_values(['g%d_minus' % col for col in xrange(A.shape[1])])
            c_plus = cpl.solution.get_values(['c%d_plus' % col for col in xrange(A.shape[1])])
            c_minus = cpl.solution.get_values(['c%d_minus' % col for col in xrange(A.shape[1])])
            
            constraint_name = 'avoid_%d' % (dimension+1)
            K[col, :] = array(c_plus) - array(c_minus)
            if LinearRegression.Rank(K, eps) < dimension+1:
                print "%d) The new mode is dependent, (dimension = %d / %d)" \
                    % (emf_counter, dimension, kernel_rank)
                K[:, dimension] *= 0.0
                continue

            print "%d) The new mode is independent, (dimension = %d / %d)" \
                % (emf_counter, dimension+1, kernel_rank)

            nonzero_values = find(K[:, dimension])
            g = min(abs(K[nonzero_values, dimension]))
            K[:, dimension] /= g
            if sum(K[:, dimension] < 0):
                K[:, dimension] *= -1.0
            cpl.linear_constraints.add(names=[constraint_name], senses='L')
            cpl.linear_constraints.set_rhs(constraint_name, len(nonzero_values)-1)
            for col in xrange(A.shape[1]):
                if g_plus[col] + g_minus[col] > 0.5:
                    cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, 1)
                    cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, 1)
                    if g_plus[col] > 0.5:
                        K[col, dimension] = c_plus[col]
                    else:
                        K[col, dimension] = -c_minus[col]
            dimension += 1
            
        return K.T
        
if __name__ == '__main__':
    A = np.matrix([[1, 0, 1, 1],[0, 1, 1, 1],[1, 1, 2, 2]])
    K = LinearRegression.FindSparseKernel(A)
    print A
    print K
    
    #y = matrix([[1],[1],[1],[1]])
    A = np.matrix([[0, 0, 0, 0, 0],[1, 0, 0, 1, 2],[2, 0, 0, 2, 4],[3,0,0,3, 6],[4,0,0,4, 8]])
    w = np.matrix([[1],[1],[1],[2],[1]])
    
    #y = A*x
    #print A
    #print x
    #w_pred, V = LinearRegression.LeastSquares(A, A*w)
    K = LinearRegression.FindSparseKernel(A)
    print A
    print K

    w_pred, V = LinearRegression.SolveLinearSystem(A, A*w)
    print w_pred
    print V
    print np.dot(A, V.T)
    

    
    x1 = np.matrix([[0,0,0,1,0]]).T
    x2 = np.matrix([[-1,0,0,-1,-2]]).T
    
    print np.linalg.norm(V*x1)<1e-10, np.dot(x1.T, w_pred)[0,0], np.dot(x1.T, w)[0,0]
    print np.linalg.norm(V*x2)<1e-10, np.dot(x2.T, w_pred)[0,0], np.dot(x2.T, w)[0,0]