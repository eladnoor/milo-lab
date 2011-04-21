#/usr/bin/python

import cvxmod
import numpy
import pylab

from cvxmod import atoms
from pygibbs.thermodynamic_constants import default_T, R


RT = R*default_T


class Pathway(object):
    """Container for doing pathway-level thermodynamic analysis."""
    
    DEFAULT_FORMATION_LB = -1e6
    DEFAULT_FORMATION_UB = 1e6
    DEFAULT_REACTION_LB = -1e3
    DEFAULT_REACTION_UB = 0.0
    
    def __init__(self, S, formation_energies):
        """Create a pathway object.
        
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            formation_energies: formation energies for the compounds
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.array format.
        """
        self.S = S
        self.formation_energies = formation_energies
        self.dG0_f = formation_energies[:, 0].tolist()
        self.Nr, self.Nc = S.shape
        
        assert len(self.dG0_f) == self.Nc

    def _RunThermoProblem(self, dgf_primes, output_var, problem):
        """Helper that runs a thermodynamic cvxmod.problem.
        
        Args:
            dgf_primes: the variable for all the transformed formation energies.
            output_var: the output variable (pCr, MTDF, etc.)
            problem: the problem object.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal output value).
        """
        problem.solve(quiet=True)   
                
        opt_val = cvxmod.value(output_var)
        opt_dgs = numpy.array(cvxmod.value(dgf_primes))        
        concentrations = pylab.exp((opt_dgs - self.formation_energies)/RT)
        return opt_dgs, concentrations, opt_val

    def _MakeFormationBounds(self, bounds=None,
                             c_range=None):
        """Make the bounds on formation energies.
        
        Args:
            bounds: a list of (ub, lb) tuples with per-compound bounds.
            c_range: the allowed range of concentrations. If not provided,
                will allow all dGf values to float between the default 
                upper and lower bounds.
        """
        # All compound activities are limited -1e6 < dGf < 1e6.
        formation_lb = [self.DEFAULT_FORMATION_LB]*self.Nc
        formation_ub = [self.DEFAULT_FORMATION_UB]*self.Nc
        
        # If a concentration range is provided, constrain
        # formation energies accordingly.
        if c_range:
            c_lower, c_upper = c_range
            make_lb = lambda x: x + RT*pylab.log(c_lower)
            make_ub = lambda x: x + RT*pylab.log(c_upper)
            formation_lb = map(make_lb, self.dG0_f)
            formation_ub = map(make_ub, self.dG0_f)

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
        """Transform all the dGf values to the given concentration.
        
        Args:
            c_mid: the concentration to transform to.
        
        Returns:
            A list of transformed dG values.
        """
        to_mid = lambda x: x + RT*pylab.log(c_mid)
        return map(to_mid, self.dG0_f)
    
    def _MakeMtdfProblem(self, c_range=(1e-6, 1e-2), bounds=None):
        """Create a cvxmod.problem for finding the Maximal Thermodynamic
        Driving Force (MTDF).
        
        Does not set the objective function... leaves that to the caller.
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A tuple (dgf_var, motive_force_var, problem_object).
        """
        # Make the formation energy variables and bounds.
        dgf_primes = cvxmod.optvar('dGf', self.Nc)
        formation_lb, formation_ub = self._MakeFormationBounds(bounds=bounds,
                                                               c_range=c_range)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr to flow forward.
        reaction_lb = cvxmod.matrix([self.DEFAULT_REACTION_LB]*self.Nr)
        reaction_ub = cvxmod.matrix([self.DEFAULT_REACTION_UB]*self.Nr)
        
        # Make the objective and problem.
        motive_force = cvxmod.optvar('B', 1)
        problem = cvxmod.problem()
        
        # Set the constraints.
        S = cvxmod.matrix(self.S)
        problem.constr.append(dgf_primes >= formation_lb)
        problem.constr.append(dgf_primes <= formation_ub)
        problem.constr.append(S*dgf_primes - motive_force <= reaction_ub)
        problem.constr.append(S*dgf_primes >= reaction_lb)
        
        return dgf_primes, motive_force, problem

    def FindMtdf(self, c_range=(1e-6, 1e-2), bounds=None):
        """Find the MTDF (Maximal Thermodynamic Driving Force).
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mtdf).
        """
        dgf_primes, motive_force, problem = self._MakeMtdfProblem(c_range, bounds)
        problem.objective = cvxmod.minimize(motive_force)
        return self._RunThermoProblem(dgf_primes, motive_force, problem)

    def FindMtdf_Regularized(self, c_range=(1e-6, 1e-2), bounds=None,
                             c_mid=1e-3,
                             max_mtdf=None):
        """Find the MTDF (Maximal Thermodynamic Driving Force).
        
        Uses l2 regularization to minimize the log difference of 
        concentrations from c_mid.
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            c_mid: the defined midpoint concentration.
            max_mtdf: the maximum value for the motive force.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mtdf).
        """
        dgf_primes, motive_force, problem = self._MakeMtdfProblem(c_range, bounds)
        dg_mids = cvxmod.matrix(self._MakeDgMids(c_mid))
        
        # Set the objective and solve.
        norm2_resid = cvxmod.norm2(dgf_primes - dg_mids)
        if max_mtdf is not None:
            problem.constr.append(motive_force <= max_mtdf)
            problem.objective = cvxmod.minimize(norm2_resid)
        else:
            problem.objective = cvxmod.minimize(motive_force + norm2_resid)

        return self._RunThermoProblem(dgf_primes, motive_force, problem)

    def _MakePcrBounds(self, c_mid, bounds):
        """Make pCr-related bounds on the formation energies.
        
        Simply: since some compounds have specific concentration
        bounds, we don't want to constrain them to be inside the pCr.
        Therefore, we need to allow their transformed formation energies
        to float in the full allowed range. 
        
        Args:
            c_mid: the defined midpoint concentration.
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A two tuple of cvxmod.matrix (lower bounds, upper bounds).
        """
        to_mid = lambda x: x + RT*pylab.log(c_mid)
        pcr_lb = self._MakeDgMids(c_mid)
        pcr_ub = self._MakeDgMids(c_mid)
        
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb:
                    pcr_lb[i] = self.DEFAULT_FORMATION_LB
                if ub:
                    pcr_ub[i] = self.DEFAULT_FORMATION_UB

        return cvxmod.matrix(pcr_lb), cvxmod.matrix(pcr_ub)

    def _MakePcrProblem(self, c_mid=1e-3, ratio=3.0,
                        bounds=None):
        """Create a cvxmod.problem for finding the pCr.
        
        Does not set the objective function... leaves that to the caller.
        
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
        reaction_lb = cvxmod.matrix([self.DEFAULT_REACTION_LB]*self.Nr)
        reaction_ub = cvxmod.matrix([self.DEFAULT_REACTION_UB]*self.Nr)
        
        # Make bounds relating pCr and formation energies.
        pcr_lb, pcr_ub = self._MakePcrBounds(c_mid, bounds)
        
        r_frac_1 = ratio / (ratio + 1.0)
        r_frac_2 = 1.0 / (ratio + 1)
        ln10 = pylab.log(10)
        
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
        return self._RunThermoProblem(dgf_primes, pc, problem)
    
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
            
        return self._RunThermoProblem(dgf_primes, pc, problem)    


if __name__ == '__main__':
    S = numpy.array([[-1,1,0],[0,-1,1]])
    cvxmod.randseed()
    dGs = numpy.array(cvxmod.randn(3,1,std=1000))
    print 'Random dGs'
    print dGs
    path = Pathway(S, dGs)
    dgs, concentrations, pcr = path.FindPcr()
    print 'Take 1', pcr
    print 'dGs'
    print dgs
    print 'concentrations'
    print concentrations
    
    dgs, concentrations, pcr = path.FindPcr_Regularized(max_pcr=pcr)
    print 'Take 2 (regularized)', pcr
    print 'dGs'
    print dgs
    print 'concentrations'
    print concentrations
