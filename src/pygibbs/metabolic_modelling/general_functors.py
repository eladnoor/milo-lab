#!/usr/bin/python

import numpy as np

from pygibbs.thermodynamic_constants import default_T, R


RT = R * default_T


class FixedVariableInjector(object):
    """Injects fixed variables into a matrix."""
    
    def __init__(self, lower_bounds, upper_bounds, initial_conds):
        """Initialize the injector.
        
        All inputs should be Numpy matrices.
        
        Args:
            lower_bounds: 1xN matrix of lower bounds on variables.
            upper_bounds: 1xN matrix of upper bounds on variables.
            initial_conds: 1xN matrix of initial conditions.
        """
        self.lb = lower_bounds
        self.ub = upper_bounds
        self.x0 = initial_conds
        self.n  = self.x0.size
        
        lb_problem_i = np.where(self.x0 < self.lb)
        ub_problem_i = np.where(self.x0 > self.ub)
        self.x0[lb_problem_i] = self.lb[lb_problem_i]
        self.x0[ub_problem_i] = self.ub[ub_problem_i]
        
        self.fixed_i = np.where(lower_bounds == upper_bounds)
        self.variable_i = np.where(lower_bounds != upper_bounds)
        self.fixed_vals = self.x0[self.fixed_i]
        
    def GetVariableLowerBounds(self):
        return self.lb[self.variable_i]
    
    def GetVariableUpperBounds(self):
        return self.ub[self.variable_i]
    
    def GetVariableInitialConds(self):
        return self.x0[self.variable_i]
    
    def __call__(self, x):
        out = np.matrix(np.zeros((1, self.n)))
        out[self.fixed_i] = self.fixed_vals
        out[self.variable_i] = x
        return out
    
    def Deject(self, injected):
        return np.array(injected[self.variable_i].flat)
    

class ConstraintChecker(object):
    """Base class for functions that check constraints.
    
    All such classes implement a __call__ method returning an np.matrix.
    """
    def __call__(self, x):
        """Returns matrix with values >= 0 where constraints are met.""" 
        raise NotImplementedError
    

class MinusDG(ConstraintChecker):
    """A functor checking that thermodynamic requirements are met."""
    
    def __init__(self, S, dG0r, injector,
                 max_dGr=0.0):
        """Initialize the MinusDG functor.
        
        Args:
            S: MxN stoichiometric matrix (Numpy matrix).
            dG0r: 1xN Numpy matrix of standard reaction energies.
            injector: FixedVariableInjector instance.
            max_dG: the maximum allowed dGr' value.
        """
        self.S = S
        self.injector = injector
        self.Ncompounds, self.Nreactions = self.S.shape
        self.dG0r = dG0r
        self.max_dGr = max_dGr
    
    def __call__(self, x):
        """Returns a matrix of -dGr so that all values are >= 0.0
           when all reactions are feasible."""
        my_x = self.injector(x)
        
        dgtag = self.dG0r + RT * my_x * self.S
        minus_max = -(dgtag - self.max_dGr)
        return minus_max
    

class BoundDiffs(ConstraintChecker):
    """Functor checking bounds."""
    
    def __init__(self, lb, ub):
        """Initialize.
        
        Args:
            lb: lower bounds.
            ub: upper bounds.
        """
        self.lb = lb
        self.ub = ub
    
    def __call__(self, x):
        """Check difference from bounds.
        
        Returns:
            A matrix where values are > 0 if bounds are met and < 0 otherwise.
        """
        lb_diff = x - self.lb
        ub_diff = self.ub - x
        return np.matrix(np.hstack([lb_diff, ub_diff]))
    

class MultiFunctionWrapper(ConstraintChecker):
    """Wraps multiple functors that return constraint-check matrices."""
    
    def __init__(self, functions):
        """Initialize the multi-functor.
        
        Args:
            functions: a list of callables.
        """
        self.functions = functions
        
    def __call__(self, x):
        """Returns the constraint check matrix 
           for all underlying functors."""
        outs = [f(x) for f in self.functions]
        out_mat = np.hstack(outs)
        
        # Required to return a 1-d ndarray.
        return np.array(out_mat.flat)