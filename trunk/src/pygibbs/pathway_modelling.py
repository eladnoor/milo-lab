#/usr/bin/python2.6

import cvxmod
from cvxmod import atoms
import math
import numpy
import pylab

from pygibbs.thermodynamic_constants import default_T, R


RT = R*default_T


class Pathway(object):
    
    def __init__(self, S, formation_energies):
        """Create a pathway object.
        
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            formation_energies: formation energies for the compounds
                in standard conditions. Should be a column vector in
                numpy.array format.
        """
        self.S = S
        self.formation_energies = formation_energies
        self.dG0_f = formation_energies[:,0].tolist()
        self.Nr, self.Nc = S.shape
        
        assert len(self.dG0_f) == self.Nc
    
    def _MakeFormationBounds(self, bounds=None):
        """Make the bounds on formation energies."""
        # All compound activities are limited -1e6 < dGf < 1e6.
        formation_lb = [-1e6]*self.Nc
        formation_ub = [1e6]*self.Nc

        # Add specific bounds for compounds with known concentrations.
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                dgf = self.dG0_f[i]
                if lb != None:
                    formation_lb[i] = dgf + RT*pylab.log(lb)
                if ub != None:
                    formation_ub[i] = dgf + RT*pylab.log(ub)
        
        return cvxmod.matrix(formation_lb), cvxmod.matrix(formation_ub)
    
    def _MakeDgMids(self, c_mid):
        to_mid = lambda x: x + RT*pylab.log(c_mid)
        return map(to_mid, self.dG0_f)
    
    def _MakePcrBounds(self, c_mid, bounds):
        """Make bounds relating to the pCr."""
        to_mid = lambda x: x + RT*pylab.log(c_mid)
        pcr_lb = self._MakeDgMids(c_mid)
        pcr_ub = self._MakeDgMids(c_mid)
        
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb:
                    pcr_lb[i] = -1e6
                if ub:
                    pcr_ub[i] = 1e6

        return cvxmod.matrix(pcr_lb), cvxmod.matrix(pcr_ub)
    
    def _MakePcrProblem(self, c_mid=1e-3, ratio=3.0,
                        bounds=None):
        """Create a cvxmod.problem for finding the pCr.
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A tuple (dgf_var, pcr_var, problem_object).
        """
        # Make the formation energy variables and bounds.
        dgf_primes = cvxmod.optvar('dGf', self.Nc)
        formation_lb, formation_ub = self._MakeFormationBounds(bounds)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr to flow forward.
        reaction_lb = cvxmod.matrix([-1e3]*self.Nr)
        reaction_ub = cvxmod.matrix([0]*self.Nr)
        
        # Make bounds relating pCr and formation energies.
        pcr_lb, pcr_ub = self._MakePcrBounds(c_mid, bounds)
        
        r_frac_1 = ratio / (ratio + 1.0)
        r_frac_2 = 1.0 / (ratio + 1)
        ln10 = math.log(10)
        
        # Make the objective and problem.
        pc = cvxmod.optvar('pC', 1)
        problem = cvxmod.problem()
        
        # Set the constraints.
        S = cvxmod.matrix(self.S)
        problem.constr.append(dgf_primes >= formation_lb)
        problem.constr.append(dgf_primes <= formation_ub)
        problem.constr.append(S*dgf_primes <= reaction_ub)
        problem.constr.append(S*dgf_primes >= reaction_lb)
        problem.constr.append(dgf_primes + r_frac_1 * RT * ln10 * pc >= pcr_lb)
        problem.constr.append(dgf_primes - r_frac_2 * RT * ln10 * pc <= pcr_ub)
        
        return dgf_primes, pc, problem
    
    def _RunPcrProblem(self, dgf_primes, pc, problem):
        """Runs the pCr problem.
        
        Args:
            dgf_primes: the variable for all the transformed formation energies.
            pc: the pCr variable.
            problem: the problem object.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal pCr value).
        """
        problem.solve(quiet=True)   
                
        opt_pc = cvxmod.value(pc)
        opt_dgs = numpy.array(cvxmod.value(dgf_primes))        
        concentrations = pylab.exp((opt_dgs - self.formation_energies)/RT)
        return opt_dgs, concentrations, opt_pc
    
    def FindPcr(self, c_mid=1e-3, ratio=3.0,
                bounds=None):
        """Compute the pCr using hip-hoptimzation!
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
 
        Returns:
            A 3 tuple (dGfs, concentrations, pCr value).
        """
        dgf_primes, pc, problem = self._MakePcrProblem(c_mid, ratio, bounds)
        
        # Set the objective and solve.
        problem.objective = cvxmod.minimize(pc)
        return self._RunPcrProblem(dgf_primes, pc, problem)
    
    def FindPcr_Regularized(self, c_mid=1e-3, ratio=3.0,
                            bounds=None, max_pcr=None):
        """Compute the pCr using hip-hoptimzation!
        
        Uses l2 regularization to minimize the log difference of 
        concentrations from c_mid.
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            max_pcr: the maximum value of the pCr to return. May be from a 
                previous, unregularized run.
 
        Returns:
            A 3 tuple (dGfs, concentrations, pCr value).
        """
        dgf_primes, pc, problem = self._MakePcrProblem(c_mid, ratio, bounds)
        dg_mids = cvxmod.matrix(self._MakeDgMids(c_mid))
        
        # Set the objective and solve.
        norm2_resid = cvxmod.norm2(dgf_primes - dg_mids)
        if max_pcr is not None:
            problem.constr.append(pc <= max_pcr)
            problem.objective = cvxmod.minimize(norm2_resid)
        else:
            problem.objective = cvxmod.minimize(pc + norm2_resid)
            
        return self._RunPcrProblem(dgf_primes, pc, problem)


if __name__ == '__main__':
    S = numpy.array([[-1,1,0],[0,-1,1]])
    cvxmod.randseed()
    dGs = numpy.array(cvxmod.randn(3,1,std=1000))
    print 'Random dGs'
    print dGs
    path = Pathway(S, dGs)
    dgs, concentrations, pcr = path.FindPcr()
    print 'Take 1', pcr
    
    dgs, concentrations, pcr = path.FindPcr_Regularized(max_pcr=pcr)
    print 'Take 2', pcr 
