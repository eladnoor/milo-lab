import cplex
from pylab import linalg, find, zeros, array, matrix
import sys

class SparseKernel(object):
    """
        Finds a sparse representation of the kernel matrix, using MILP
        to iterate the Fundamental Modes of the matrix.
    
        Input:
            a (n x m) matrix A, whose rank is r.
        Return:
            a (m x m-r) matrix K that will span the kernel of A, i.e.:
            span(K) = {x | Ax = 0}
    """

    class LinearProgrammingException(Exception):
        pass

    @staticmethod
    def Rank(A, eps=1e-10):
        _U, s, _V = linalg.svd(A, full_matrices=False)
        return len(find(s > eps))
    
    def __init__(self, A):
        self.upper_bound = 1000
        self.eps = 1e-10
        self.A = A
        self.cpl = cplex.Cplex()
        self.cpl.set_problem_name('find_kernel')
        self.cpl.set_log_stream(None)
        self.cpl.set_results_stream(None)
        self.cpl.set_warning_stream(None)
        self.kernel_rank = self.A.shape[1] - SparseKernel.Rank(self.A, self.eps)
        self.CreateAllVariables()
        self.AddLinearConstraints()

    def CreateAllVariables(self):
        """ 
            create 4 variables for each column: 
            positive & negative real values and positive & negative indicators
        """
        for col in xrange(self.A.shape[1]):
            self.cpl.variables.add(names=['c%d_plus' % col], 
                                   lb=[0], ub=[self.upper_bound])
            self.cpl.variables.add(names=['g%d_plus' % col], types='B')
    
            # c_plus - M*g_plus <= 0
            constraint_name = 'g%d_plus_bound_upper' % col
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_plus' % col, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, -self.upper_bound)
    
            # c_plus - g_plus >= 0
            constraint_name = 'g%d_plus_bound_lower' % col
            self.cpl.linear_constraints.add(names=[constraint_name], senses='G', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_plus' % col, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, -1)
    
            self.cpl.variables.add(names=['c%d_minus' % col], lb=[0], ub=[self.upper_bound])
            self.cpl.variables.add(names=['g%d_minus' % col], types='B')
    
            # c_minus - M*g_minus <= 0
            constraint_name = 'g%d_minus_bound_upper' % col
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_minus' % col, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, -self.upper_bound)
    
            # c_minus - g_minus >= 0
            constraint_name = 'g%d_minus_bound_lower' % col
            self.cpl.linear_constraints.add(names=[constraint_name], senses='G', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'c%d_minus' % col, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, -1)
            
            # g_plus + g_minus <= 1
            constraint_name = 'g%d_bound' % col
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[1])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, 1)


    def AddLinearConstraints(self):
        # add the linear constraints corresponding to the values in the matrix A
        for row in xrange(self.A.shape[0]):
            constraint_name = 'r%d' % row
            self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
            for col in xrange(self.A.shape[1]):
                if self.A[row,col] != 0:
                    self.cpl.linear_constraints.set_coefficients(constraint_name, 
                        'c%d_plus' % col, self.A[row,col])
                    self.cpl.linear_constraints.set_coefficients(constraint_name,
                        'c%d_minus' % col, -self.A[row,col])
    
        all_gammas = ['g%d_plus' % col for col in xrange(self.A.shape[1])] + \
                     ['g%d_minus' % col for col in xrange(self.A.shape[1])]
        
        # Set the objective function: minimizing the sum of gammas
        self.cpl.objective.set_linear([(g, 1) for g in all_gammas])
    
        # Exclude the trivial solution (all-zeros)
        self.cpl.linear_constraints.add(names=['avoid_0'], senses='G', rhs=[1])
        for g in all_gammas:
            self.cpl.linear_constraints.set_coefficients('avoid_0', g, 1)

    def GetSolution(self):
            try:
                self.cpl.solve()
            except cplex.exceptions.CplexSolverError as e:
                raise SparseKernel.LinearProgrammingException(str(e))
            
            sol = self.cpl.solution
            if self.cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.MIP_optimal:
                raise SparseKernel.LinearProgrammingException("No more EMFs")
            
            g_plus = array(sol.get_values(['g%d_plus' % col for col in xrange(self.A.shape[1])]))
            g_minus = array(sol.get_values(['g%d_minus' % col for col in xrange(self.A.shape[1])]))
            
            coeffs = array(sol.get_values(['c%d_plus' % col for col in xrange(self.A.shape[1])])) - \
                     array(sol.get_values(['c%d_minus' % col for col in xrange(self.A.shape[1])]))
            
            return g_plus, g_minus, coeffs
            
    def ExcludeSolutionVector(self, g_plus, g_minus, constraint_name):
        self.cpl.linear_constraints.add(names=[constraint_name], senses='L')
        nonzero_counter = 0
        for col in xrange(self.A.shape[1]):
            if g_plus[col] > 0.5:
                self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_plus' % col, 1)
                nonzero_counter += 1
            elif g_minus[col] > 0.5:
                self.cpl.linear_constraints.set_coefficients(constraint_name, 'g%d_minus' % col, 1)
                nonzero_counter += 1
        
        self.cpl.linear_constraints.set_rhs(constraint_name, nonzero_counter-1)

    def __len__(self):
        return self.kernel_rank

    def __iter__(self):
        self.dimension = 0
        self.emf_counter = 0
        self.K = zeros((self.A.shape[1], len(self)))
        return self
    
    def next(self):
        if self.dimension == self.K.shape[1]:
            raise StopIteration
        
        while True:
            self.emf_counter += 1
            
            g_plus, g_minus, coeffs = self.GetSolution()
            self.ExcludeSolutionVector(g_plus, g_minus, 'avoid_%d_plus' % self.emf_counter)
            self.ExcludeSolutionVector(g_minus, g_plus, 'avoid_%d_minus' % self.emf_counter)
            
            nonzero_indices = find(g_plus > 0.5).tolist() + find(g_minus > 0.5).tolist()
            self.K[nonzero_indices, self.dimension] = coeffs[nonzero_indices]
            
            if SparseKernel.Rank(self.K, self.eps) < self.dimension+1:
                self.K[:, self.dimension] = 0.0
            else:
                # normalize the kernel vector so that it will have nice coefficients
                nonzero_values = find(self.K[:, self.dimension])
                g = min(abs(self.K[nonzero_values, self.dimension]))
                self.K[:, self.dimension] /= g
                if sum(self.K[:, self.dimension] < 0.0):
                    self.K[:, self.dimension] *= -1.0
                
                v = self.K[:, self.dimension]
                self.dimension += 1
                return v
        
    def Solve(self):
        for _ in self:
            pass
        return self.K

if __name__ == '__main__':
    A = matrix([[1, 0, 1, 1],[0, 1, 1, 1],[1, 1, 2, 2]])
    K = SparseKernel(A).Solve()
    print A
    print K
    print A*K
