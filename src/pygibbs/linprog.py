import cvxopt
import cplex
import pylab

def linprog_glpk(f, A, b, lb=[], ub=[]):
    """
        The constraints are:
            A*x <= b
            lb <= x <= ub
        
        The optimization rule is:
            minimize f*x
        
        All other parameters (f, b, lb, ub) should be column vectors, with the following sizes:
            A, Aeq - Nr x Nc matrix
            f      - Nc x 1 matrix
            b, beq - Nr x 1 matrix
            lb     - list of pairs (column_index, lower_bound)
            ub     - list of pairs (column_index, upper_bound)
            
        Returns:
            x - the solution to the linear problem
        or
            None - if there is no solution
    """
    (Nr, Nc) = A.shape
    if (f.shape[0] != Nc or f.shape[1] != 1):
        raise Exception("linprog: 'f' must be a column vector whose length matches the number of columns in 'A'")
    if (b.shape[0] != Nr or b.shape[1] != 1):
        raise Exception("linprog: 'b' must be a column vector whose length matches the number of rows in 'A'")
    
    for (c, bound) in lb:
        row_lower = pylab.zeros(Nc)
        row_lower[c] = -1
        A = pylab.vstack([A, row_lower])
        b = pylab.vstack([b, [-bound]])
    for (c, bound) in ub:
        row_upper = pylab.zeros(Nc)
        row_upper[c] = 1
        A = pylab.vstack([A, row_upper])
        b = pylab.vstack([b, [bound]])
    
    # make sure the values are floating point (not integer)
    # and convert the matrices to cvxopt format
    c = cvxopt.matrix(f*1.0)
    G = cvxopt.matrix(A*1.0)
    h = cvxopt.matrix(b*1.0)
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['LPX_K_MSGLEV'] = 0
    sol = cvxopt.solvers.lp(c, G, h, solver='glpk')
    if (not sol['x']):
        return None
    else:
        return pylab.matrix(sol['x'])


def linprog(f, A, b, A_eq=None, b_eq=None, lb=None, ub=None, log_stream=None):
    """
        The constraints are:
            A*x <= b
            Aeq*x = beq  - optional
            lb <= x <= ub
        
        The optimization rule is:
            minimize f*x
        
        All other parameters (f, b, lb, ub) should be column vectors, with the following sizes:
            f   - Nc x 1 matrix
            A   - Nr x Nc matrix
            b   - Nr x 1 matrix
            Aeq - Nr_eq x Nc matrix
            beq - Nr_eq x 1 matrix
            lb  - Nc x 1 matrix - or - list of pairs (column_index, lower_bound)
            ub  - Nc x 1 matrix - or - list of pairs (column_index, upper_bound)
            
        Returns:
            x - the solution to the linear problem (Nc x 1 matrix)
        or
            None - if there is no solution
    """
    if (f.shape[1] != 1):
        raise Exception("linprog: 'f' must be a column vector")
    Nc = f.shape[0]
    
    if (A != None):
        if (A.shape[1] != Nc):
            raise Exception("'A' must have the same number of columns as 'f'")
        Nr = A.shape[0]
    else:
        Nr = 0

    if (A_eq != None):
        if (A_eq.shape[1] != Nc):
            raise Exception("'Aeq' must have the same number of columns as 'f'") 
        Nr_eq = A_eq.shape[0]
    else:
        Nr_eq = 0

    if (Nr == 0 and Nr_eq == 0):
        raise Exception("linprog: neither A nor Aeq were provided")
    
    if (b != None and b.shape != (Nr, 1)):
        raise Exception("linprog: 'b' must be a column vector whose length matches the number of rows in 'A'")
    if (b_eq != None and b_eq.shape != (Nr_eq, 1)):
        raise Exception("linprog: 'b_eq' must be a column vector whose length matches the number of rows in 'A_eq'")

    if (lb != None and type(lb) == type(pylab.array([]))):
        if (lb.shape != (Nc, 1)):
            raise Exception("linprog: 'lb' must be a column vector whose length matches the number of columns in 'A'")
        lb = [(c, lb[c,0]) for c in range(Nc)]

    if (ub != None and type(ub) == type(pylab.array([]))):
        if (ub.shape != (Nc, 1)):
            raise Exception("linprog: 'lb' must be a column vector whose length matches the number of columns in 'A'")
        ub = [(c, ub[c,0]) for c in range(Nc)]

    cpl = cplex.Cplex()
    if (log_stream != None):
        cpl.set_log_stream(log_stream)
    else:
        cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    #cpl.set_warning_stream(None)
    
    cpl.set_problem_name('LP')
    cpl.variables.add(names=["c%d" % c for c in range(Nc)])
    if (lb != None):
        cpl.variables.set_lower_bounds(lb)
    if (ub != None):
        cpl.variables.set_upper_bounds(ub)

    cpl.objective.set_linear([(c, f[c,0]) for c in range(Nc)])
    
    if (Nr > 0):
        for r in range(Nr):
            cpl.linear_constraints.add(senses='L', names=["r%d" % r])
            for c in pylab.find(A[r,:] != 0):
                cpl.linear_constraints.set_coefficients(r, c, A[r,c])
        cpl.linear_constraints.set_rhs([(r, b[r,0]) for r in range(Nr)])

    if (Nr_eq > 0):
        for r in range(Nr_eq):
            cpl.linear_constraints.add(senses='E', names=["r%d" % (Nr+r)])
            for c in pylab.find(A_eq[r,:] != 0):
                cpl.linear_constraints.set_coefficients(Nr+r, c, A_eq[r,c])
        cpl.linear_constraints.set_rhs([(Nr+r, b_eq[r,0]) for r in range(Nr_eq)])

    cpl.solve()
    if (cpl.solution.get_status() != cplex.callbacks.SolveCallback.status.optimal):
        return None
    else:
        return pylab.matrix(cpl.solution.get_values()).T

if (__name__ == '__main__'):
    f = pylab.matrix([1, 1]).T
    A = pylab.matrix([[-1, -2],[-2, -1],[0, -1]])
    b = pylab.matrix([-3,-5,0]).T
    A_eq = pylab.matrix([[-1, -2],[-2, -1]])
    b_eq = pylab.matrix([-3,-5]).T
    print linprog(f, A, b, None, None)
    print linprog(f, None, None, A_eq, b_eq)
