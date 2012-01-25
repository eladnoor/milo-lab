#/usr/bin/python

import cvxmod
import numpy as np
import matplotlib.pyplot as plt

from pygibbs.thermodynamic_constants import default_T, R
from matplotlib.font_manager import FontProperties
from pygibbs.kegg import Kegg
from pygibbs.kegg_reaction import Reaction
from pygibbs.metabolic_modelling.mtdf_optimizer import UnsolvableConvexProblemException
import types
import sys

RT = R * default_T

class DeltaGNormalization:
    TIMES_FLUX = 1 # motivation is having the most uniform entropy production
    DIVIDE_BY_FLUX = 2 # motivation is requiring more force in reaction that have more flux
    SIGN_FLUX = 3 # motivation is putting limits on the allowed backward/forward fluxes
    
    DEFAULT = SIGN_FLUX

class Pathway(object):
    """Container for doing pathway-level thermodynamic analysis."""
    
    DEFAULT_FORMATION_LB = -1e6
    DEFAULT_FORMATION_UB = 1e6
    DEFAULT_REACTION_LB = -1e3
    DEFAULT_REACTION_UB = 0.0
    DEFAULT_C_RANGE = (1e-6, 0.1)
    DEFAULT_PHYSIOLOGICAL_CONC = 1e-3
    
    def __init__(self, S,
                 formation_energies=None, reaction_energies=None, fluxes=None):
        """Create a pathway object.
        
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            formation_energies: the Gibbs formation energy for the compounds
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.array format.
            reaction_energies: the change in Gibbs energy for the reactions
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.array format.
            fluxes: the list of relative fluxes through each of the reactions.
                By default, all fluxes are 1.
        """
        if formation_energies is None and reaction_energies is None:
            raise ValueError("In order to use 'Pathway' xeither formation xor "
                             "reaction energies must be provided.")
        if formation_energies is not None and reaction_energies is not None:
            raise ValueError("In order to use 'Pathway' xeither formation xor "
                             "reaction energies must be provided.")
        
        self.S = S
        self.dG0_f_prime = formation_energies
        self.dG0_r_prime = reaction_energies
        self.Nr, self.Nc = S.shape
        
        if self.dG0_f_prime is not None:
            assert self.dG0_f_prime.shape[0] == self.Nc
            self.dG0_r_prime = self.CalculateReactionEnergies(self.dG0_f_prime)
        else:
            assert self.dG0_r_prime.shape[0] == self.Nr

        if fluxes is None:
            self.fluxes = [1] * self.Nr
        else:
            assert len(fluxes) == self.Nr
            self.fluxes = fluxes

    def CalculateReactionEnergies(self, dG_f):
        if np.isnan(dG_f).any():
            # if there are NaN values in dG_f, multiplying the matrices will not
            # work, since NumPy will not convert 0*NaN into 0 in the sum. Therefore,
            # the multiplication must be done explicitly and using only the nonzero
            # stoichiometric coefficients and their corresponding dG_f. 
            dG_r = np.zeros((self.Nr, 1))
            for r in xrange(self.Nr):
                reactants = np.nonzero(self.S[r, :])[0]
                dG_r[r, 0] = np.dot(self.S[r, reactants], dG_f[reactants])
            return dG_r
        else:
            return np.dot(self.S, dG_f)

    def CalculateReactionEnergiesUsingConcentrations(self, concentrations):
        log_conc = np.log(concentrations)
        if np.isnan(self.dG0_r_prime).any(): # see CalculateReactionEnergies
            dG_r_prime = np.zeros((self.Nr, 1))
            for r in xrange(self.Nr):
                reactants = np.nonzero(self.S[r, :])[0]
                dG_r_prime[r, 0] = self.dG0_r_prime[r, 0] + \
                    RT * np.dot(self.S[r, reactants], log_conc[reactants])
            return dG_r_prime
        else:
            return self.dG0_r_prime + RT * np.dot(self.S, log_conc)

    def GetPhysiologicalConcentrations(self, bounds=None):
        conc = np.ones((self.Nc, 1)) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None and ub is not None:
                    if not (lb < conc[i, 0] < ub):
                        conc[i, 0] = np.sqrt(lb * ub)
        
        return conc

    def _RunThermoProblem(self, ln_conc, output_var, problem):
        """Helper that runs a thermodynamic cvxmod.problem.
        
        Args:
            dgf_primes: the variable for all the transformed formation energies.
            output_var: the output variable (pCr, MTDF, etc.)
            problem: the problem object.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal output value).
        """
        try:
            status = problem.solve(quiet=True)
        except ValueError as e:
            sys.stderr.write(str(problem))
            raise UnsolvableConvexProblemException(str(e), problem)
        if status != 'optimal':
            raise UnsolvableConvexProblemException(status, problem)
            
        opt_val = cvxmod.value(output_var)
        opt_ln_conc = np.array(cvxmod.value(ln_conc))        
        concentrations = np.exp(opt_ln_conc)
        return opt_ln_conc, concentrations, opt_val

    def _MakeLnConcentratonBounds(self, bounds=None, c_range=None):
        """Make bounds on logarithmic concentrations."""
        _c_range = c_range or self.DEFAULT_C_RANGE
        c_lower, c_upper = c_range
        ln_conc_lb = cvxmod.matrix([np.log(c_lower)] * self.Nc, (self.Nc, 1))
        ln_conc_ub = cvxmod.matrix([np.log(c_upper)] * self.Nc, (self.Nc, 1))
        
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None:
                    ln_conc_lb[i, 0] = np.log(lb)
                if ub is not None:
                    ln_conc_ub[i, 0] = np.log(ub)
        
        return ln_conc_lb, ln_conc_ub
            
    def _MakeFormationBounds(self, bounds=None, c_range=None):
        """Make the bounds on formation energies.
        
        Args:
            bounds: a list of (ub, lb) tuples with per-compound bounds.
            c_range: the allowed range of concentrations. If not provided,
                will allow all dGf values to float between the default 
                upper and lower bounds.
        """
        assert self.dG0_f_prime is not None
        
        if c_range:
            # If a concentration range is provided, constrain
            # formation energies accordingly.
            c_lower, c_upper = c_range
            formation_lb = cvxmod.matrix(self.dG0_f_prime + RT * np.log(c_lower))
            formation_ub = cvxmod.matrix(self.dG0_f_prime + RT * np.log(c_upper))
        else:
            # Otherwise, all compound activities are limited -1e6 < dGf < 1e6.
            formation_lb = cvxmod.matrix(self.DEFAULT_FORMATION_LB, size=(self.Nc, 1))
            formation_ub = cvxmod.matrix(self.DEFAULT_FORMATION_UB, size=(self.Nc, 1))

        # Override the specific for compounds with known concentrations.
        if bounds:
            for i, bound in enumerate(bounds):
                dgf = self.dG0_f_prime[i, 0]
                if np.isnan(dgf):
                    continue
                lb, ub = bound
                if lb is not None:
                    formation_lb[i, 0] = dgf + RT * np.log(lb)
                if ub is not None:
                    formation_ub[i, 0] = dgf + RT * np.log(ub)

        lb_nan_indices = np.isnan(formation_lb).nonzero()[0].tolist()
        ub_nan_indices = np.isnan(formation_ub).nonzero()[0].tolist()
        formation_lb[lb_nan_indices] = self.DEFAULT_FORMATION_LB
        formation_ub[ub_nan_indices] = self.DEFAULT_FORMATION_UB
        
        return formation_lb, formation_ub
    
    def _MakeDgMids(self, c_mid):
        """Transform all the dGf values to the given concentration.
        
        Args:
            c_mid: the concentration to transform to.
        
        Returns:
            A list of transformed dG values.
        """
        assert self.dG0_f_prime is not None

        to_mid = lambda x: x + RT * np.log(c_mid)
        return map(to_mid, self.dG0_f_prime[:, 0].tolist())

    def _GetTotalReactionEnergy(self, c_range=(1e-6, 1e-2), bounds=None):
        """
            Calculates the overall reaction and finds the minimal and maximal
            values for its deltaG', for the given concentration bounds.
        """
        problem = cvxmod.problem()

        # Define and apply the constraints on the concentrations
        ln_conc = cvxmod.optvar('lnC', self.Nc, 1)
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds(bounds=bounds,
                                                                c_range=c_range)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)
        
        # find the row vector describing the overall stoichiometry
        S = cvxmod.matrix(self.S)
        f = cvxmod.matrix(np.array(self.fluxes)).T
        g0 = cvxmod.matrix(self.dG0_r_prime)
        g = g0 + RT * S * ln_conc
        total_g = f * g
        
        problem.objective = cvxmod.maximize(total_g)
        _, _, max_total_g = self._RunThermoProblem(ln_conc, total_g, problem)

        problem.objective = cvxmod.minimize(total_g)
        _, _, min_total_g = self._RunThermoProblem(ln_conc, total_g, problem)
        
        return min_total_g, max_total_g
    
    def _GetMaxReactionEnergy(self, c_range=(1e-6, 1e-2), bounds=None, min_driving_force=0):
        problem = cvxmod.problem()

        # Define and apply the constraints on the concentrations
        ln_conc = cvxmod.optvar('lnC', self.Nc, 1)
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds(bounds=bounds,
                                                                c_range=c_range)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)
        
        # find the row vector describing the overall stoichiometry
        S = cvxmod.matrix(self.S)
        f = cvxmod.matrix(np.array(self.fluxes)).T
        g0 = cvxmod.matrix(self.dG0_r_prime)
        g = g0 + RT * S * ln_conc
        problem.constr.append(g <= -min_driving_force)

        total_g = f * g
        
        problem.objective = cvxmod.maximize(total_g)
        return self._RunThermoProblem(ln_conc, total_g, problem)
            
    def _MakeMtdfProblem(self, c_range=(1e-6, 1e-2), bounds=None,
                         normalization=DeltaGNormalization.DEFAULT):
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
        problem = cvxmod.problem()

        # Define and apply the constraints on the concentrations
        ln_conc = cvxmod.optvar('lnC', self.Nc)
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds(bounds=bounds,
                                                                c_range=c_range)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)

        # Make the objective
        motive_force_lb = cvxmod.optvar('B', 1)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium. 
        S = cvxmod.matrix(self.S)
        dg0r_primes = cvxmod.matrix(self.dG0_r_prime)
        for i, flux in enumerate(self.fluxes):
            
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(self.dG0_r_prime[i, 0]):
                continue
            
            curr_dgr = dg0r_primes[i, 0] + RT * S[i, :] * ln_conc
            if flux == 0:
                problem.constr.append(curr_dgr == 0)
            else:
                if normalization == DeltaGNormalization.DIVIDE_BY_FLUX:
                    motive_force = -curr_dgr * (1.0 / flux)
                elif normalization == DeltaGNormalization.TIMES_FLUX:
                    motive_force = -curr_dgr * flux
                elif normalization == DeltaGNormalization.SIGN_FLUX:
                    motive_force = -curr_dgr * np.sign(flux)
                else:
                    raise ValueError("bad value for normalization method: "
                                     + str(normalization))

                problem.constr.append(motive_force >= motive_force_lb)
                #problem.constr.append(motive_force <= -self.DEFAULT_REACTION_LB)
        
        return ln_conc, motive_force_lb, problem

    def _FindMtdf(self, c_range=(1e-6, 1e-2), bounds=None,
                  normalization=DeltaGNormalization.DEFAULT):
        """Find the MTDF (Maximal Thermodynamic Driving Force).
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mtdf).
        """
        ln_conc, motive_force_lb, problem = self._MakeMtdfProblem(
                                                c_range, bounds, normalization)
        problem.objective = cvxmod.maximize(motive_force_lb)
        return self._RunThermoProblem(ln_conc, motive_force_lb, problem)

    def FindMtdf_Regularized(self, c_range=(1e-6, 1e-2), bounds=None,
                             c_mid=1e-3,
                             min_mtdf=None,
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
        if max_mtdf is not None and min_mtdf is not None:
            problem.constr.append(motive_force <= max_mtdf)
            problem.constr.append(motive_force >= min_mtdf)
            problem.objective = cvxmod.minimize(norm2_resid)
        elif max_mtdf is not None:
            problem.constr.append(motive_force <= max_mtdf)
            problem.objective = cvxmod.minimize(norm2_resid)
        elif min_mtdf is not None:
            problem.constr.append(motive_force >= min_mtdf)
            problem.objective = cvxmod.minimize(norm2_resid)
        else:
            problem.objective = cvxmod.minimize(motive_force + norm2_resid)

        return self._RunThermoProblem(dgf_primes, motive_force, problem)

    def FindMTDF_OptimizeConcentrations(self, c_range=1e-3,
                                        bounds=None, c_mid=1e-3):
        """Optimize concentrations at optimal pCr.
        
        Runs two rounds of optimization to find "optimal" concentrations
        at the optimal MTDF. First finds the globally optimal MTDF.
        Then minimizes the l2 norm of deviations of log concentrations
        from c_mid given the optimal MTDF.

        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            c_mid: the median concentration.
 
        Returns:
            A 3 tuple (dGfs, concentrations, MTDF value).
        """
        _, _, opt_mtdf = self._FindMtdf(c_range, bounds)
        return self.FindMtdf_Regularized(c_range, bounds, c_mid,
                                         max_mtdf=opt_mtdf)

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
                        bounds=None, max_reaction_dg=None):
        """Create a cvxmod.problem for finding the pCr.
        
        Does not set the objective function... leaves that to the caller.
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            max_reaction_dg: the maximum reaction dG allowed. Can be used to 
                enforce more motive force than feasibility.
                
        Returns:
            A tuple (dgf_var, pcr_var, problem_object).
        """
        # Make the formation energy variables and bounds.
        dgf_primes = cvxmod.optvar('dGf', self.Nc)
        formation_lb, formation_ub = self._MakeFormationBounds(bounds)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr to flow forward.
        r_ub = max_reaction_dg or self.DEFAULT_REACTION_UB
        reaction_lb = cvxmod.matrix([self.DEFAULT_REACTION_LB] * self.Nr)
        reaction_ub = cvxmod.matrix([r_ub] * self.Nr)
        
        # Make bounds relating pCr and formation energies.
        pcr_lb, pcr_ub = self._MakePcrBounds(c_mid, bounds)
        
        r_frac_1 = ratio / (ratio + 1.0)
        r_frac_2 = 1.0 / (ratio + 1)
        ln10 = np.log(10)
        
        # Make the objective and problem.
        pc = cvxmod.optvar('pC', 1)
        problem = cvxmod.problem()
        
        # Set the constraints.
        S = cvxmod.matrix(self.S)
        problem.constr.append(dgf_primes >= formation_lb)
        problem.constr.append(dgf_primes <= formation_ub)
        problem.constr.append(S * dgf_primes <= reaction_ub)
        problem.constr.append(S * dgf_primes >= reaction_lb)
        problem.constr.append(dgf_primes + r_frac_1 * RT * ln10 * pc >= pcr_lb)
        problem.constr.append(dgf_primes - r_frac_2 * RT * ln10 * pc <= pcr_ub)
        
        return dgf_primes, pc, problem
    
    def FindPcr(self, c_mid=1e-3, ratio=3.0,
                bounds=None, max_reaction_dg=None):
        """Compute the pCr using hip-hoptimzation!
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            max_reaction_dg: the maximum reaction dG allowed. Can be used to 
                enforce more motive force than feasibility.
 
        Returns:
            A 3 tuple (dGfs, concentrations, pCr value).
        """
        dgf_primes, pc, problem = self._MakePcrProblem(c_mid, ratio, bounds,
                                                       max_reaction_dg)
        
        # Set the objective and solve.
        problem.objective = cvxmod.minimize(pc)
        return self._RunThermoProblem(dgf_primes, pc, problem)
    
    def FindPcr_Regularized(self, c_mid=1e-3, ratio=3.0,
                            bounds=None, max_reaction_dg=None, max_pcr=None):
        """Compute the pCr using hip-hoptimzation!
        
        Uses l2 regularization to minimize the log difference of 
        concentrations from c_mid.
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            max_reaction_dg: the maximum reaction dG allowed. Can be used to 
                enforce more motive force than feasibility.
            max_pcr: the maximum value of the pCr to return. May be from a 
                previous, unregularized run.
 
        Returns:
            A 3 tuple (dGfs, concentrations, pCr value).
        """
        dgf_primes, pc, problem = self._MakePcrProblem(c_mid, ratio, bounds,
                                                       max_reaction_dg)
        dg_mids = cvxmod.matrix(self._MakeDgMids(c_mid))
        
        # Set the objective and solve.
        norm2_resid = cvxmod.norm2(dgf_primes - dg_mids)
        if max_pcr is not None:
            problem.constr.append(pc <= max_pcr)
            problem.objective = cvxmod.minimize(norm2_resid)
        else:
            problem.objective = cvxmod.minimize(pc + norm2_resid)
            
        return self._RunThermoProblem(dgf_primes, pc, problem)    

    def FindPcr_OptimizeConcentrations(self, c_mid=1e-3, ratio=3.0,
                                       bounds=None, max_reaction_dg=None):
        """Optimize concentrations at optimal pCr.
        
        Runs two rounds of optimization to find "optimal" concentrations
        at the optimal pCr. First finds the globally optimal pCr.
        Then minimizes the l2 norm of deviations of log concentrations
        from c_mid given the optimal pCr.

        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            max_reaction_dg: the maximum reaction dG allowed. Can be used to 
                enforce more motive force than feasibility.
 
        Returns:
            A 3 tuple (dGfs, concentrations, pCr value).
        """
        _, _, opt_pcr = self.FindPcr(c_mid, ratio, bounds, max_reaction_dg)
        return self.FindPcr_Regularized(c_mid, ratio, bounds,
                                        max_reaction_dg=max_reaction_dg,
                                        max_pcr=opt_pcr)

    def FindPcrEnzymeCost(self, c_mid=1e-3, ratio=3.0,
                          bounds=None, max_reaction_dg=None,
                          fluxes=None, km=1e-5, kcat=1e2):
        """Find the enzymatic cost at the pCr.
        
        TODO(flamholz): account for enzyme mass or # of amino acids.
        
        Assumption: all enzymes have the same kinetics for all substrates.
        Lowest substrate concentration is limiting. All other substrates 
        are ignorable. Enzymes follow Michaelis-Menten kinetics w.r.t. limiting
        substrate.
        
        Args:
            c_mid: the median concentration.
            ratio: the ratio between the distance of the upper bound
                from c_mid and the lower bound from c_mid (in logarithmic scale).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            fluxes: the relative flux through each reaction in the pathway.
            cofactors: a mapping from compound index to boolean indicating
                whether it is a cofactor.
            km: the Km value to use for all enzymes.
            kcat: the kcat value to use for all enzymes.
        
        Returns:
            A 3-tuple (total cost, per enzyme cost, concentrations).
        """
        _dgs, concentrations, _opt_pcr = self.FindPcr_OptimizeConcentrations(
            c_mid, ratio, bounds, max_reaction_dg)
        
        fluxes = fluxes or [1.0] * self.Nr
        costs = []
        for i in xrange(self.Nr):
            rxn = self.S[i, :]
            substrate_indices = np.where(rxn < 0)[0]
            substrate_concentrations = concentrations[substrate_indices, 0]
            
            # Min concentration is rate-limiting.
            s = min(substrate_concentrations)
            v = fluxes[i]
            enzyme_units = v * (km + s) / (kcat * s)
            costs.append(enzyme_units)
        
        return sum(costs), costs, concentrations

    def _MakeMinimalFeasbileConcentrationProblem(self, bounds=None, c_range=(1e-6, 1e-2)):
        problem = cvxmod.problem()

        # Define and apply the constraints on the concentrations
        ln_conc = cvxmod.optvar('lnC', self.Nc)
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds(bounds=bounds,
                                                                c_range=c_range)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)

        # find the row vector describing the overall stoichiometry
        S = cvxmod.matrix(self.S)
        dg0r_primes = cvxmod.matrix(self.dG0_r_prime)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium. 
        for i, flux in enumerate(self.fluxes):
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(self.dG0_r_prime[i, 0]):
                continue
            
            curr_dgr = dg0r_primes[i, 0] + RT * S[i, :] * ln_conc

            if flux != 0:
                problem.constr.append(curr_dgr * np.sign(flux) <= self.DEFAULT_REACTION_UB)
                problem.constr.append(curr_dgr * np.sign(flux) >= self.DEFAULT_REACTION_LB)
            else:
                problem.constr.append(curr_dgr == 0)
        
        # Set the constraints.
        return ln_conc, problem

    def FindMinimalFeasibleConcentration(self, index_to_minimize,
                                         bounds=None, c_range=(1e-6, 1e-2)):
        """
            Compute the smallest ratio between two concentrations which makes the pathway feasible.
            All other compounds except these two are constrained by 'bounds' or unconstrained at all.
        
            Arguments:
                index_to_minimize - the column index of the compound whose concentration 
                                    is to be minimized
        
            Returns:
                dGs, concentrations, target-concentration
        """
        ln_conc, problem = self._MakeMinimalFeasbileConcentrationProblem(bounds, c_range)
        problem.objective = cvxmod.minimize(ln_conc[index_to_minimize]) 
        return self._RunThermoProblem(ln_conc, ln_conc[index_to_minimize], problem)

    def _MakeMinimumFeasbileConcentrationsProblem(self, bounds=None,
                                                  c_range=(1e-6, 1e-2)):
        """Creates the cvxmod.problem for finding minimum total concentrations.
        
        Returns:
            Two tuple (ln_concentrations var, problem).
        """
        assert self.dG0_f_prime is not None
        
        ln_concentrations = cvxmod.optvar('ln_concs', self.Nc)
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds(bounds=bounds,
                                                                c_range=c_range)
        
        # Make the objective and problem.
        problem = cvxmod.problem()
        S = cvxmod.matrix(self.S)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium.
        dgf_primes = (RT * ln_concentrations) + cvxmod.matrix(self.dG0_f_prime)
        for i, flux in enumerate(self.fluxes):
            if flux > 0:
                problem.constr.append(S[i, :] * dgf_primes <= self.DEFAULT_REACTION_UB)
                problem.constr.append(S[i, :] * dgf_primes >= self.DEFAULT_REACTION_LB)
            elif flux == 0:
                problem.constr.append(S[i, :] * dgf_primes == 0)
            else:
                problem.constr.append(S[i, :] * dgf_primes >= -self.DEFAULT_REACTION_UB)
                problem.constr.append(S[i, :] * dgf_primes <= -self.DEFAULT_REACTION_LB)
        
        # Set the constraints.
        problem.constr.append(ln_concentrations >= ln_conc_lb)
        problem.constr.append(ln_concentrations <= ln_conc_ub)
        return ln_concentrations, problem

    def FindMinimumFeasibleConcentrations(self, bounds=None):
        """Use the power of convex optimization!
        
        minimize sum (concentrations)
        
        we can do this by using ln(concentration) as variables and leveraging 
        the convexity of exponentials. 
        
        min sum (exp(ln(concentrations)))
        """
        assert self.dG0_f_prime is not None
        
        ln_concentrations, problem = self._MakeMinimumFeasbileConcentrationsProblem(bounds=bounds)
        problem.objective = cvxmod.minimize(cvxmod.sum(cvxmod.atoms.exp(ln_concentrations)))
        
        status = problem.solve(quiet=True)
        if status != 'optimal':
            raise UnsolvableConvexProblemException(status, problem)
            
        opt_val = cvxmod.value(problem.objective)
        ln_concentrations = np.array(cvxmod.value(ln_concentrations))
        concentrations = np.exp(ln_concentrations)
        opt_dgs = RT * ln_concentrations + self.dG0_f_prime
        return opt_dgs, concentrations, opt_val

    def FindKineticOptimum(self, fluxes=None, km=1e-5, kcat=1e2):
        """Use the power of convex optimization!
        
        minimize sum (protein cost)
        
        we can do this by using ln(concentration) as variables and leveraging 
        the convexity of exponentials.         
        """
        assert self.dG0_f_prime is not None
        
        ln_concentrations, problem = self._MakeMinimumFeasbileConcentrationsProblem()
        #scaled_fluxes = cvxmod.matrix(fluxes or [1.0]*self.Nr) * (km/kcat)
        problem.objective = cvxmod.minimize(cvxmod.sum(cvxmod.atoms.exp(-ln_concentrations)))
        
        print problem.classify()
        status = problem.solve(quiet=False)
        if status != 'optimal':
            
            raise UnsolvableConvexProblemException(status, problem)
            
        opt_val = cvxmod.value(problem.objective)
        ln_concentrations = np.array(cvxmod.value(ln_concentrations))
        concentrations = np.exp(ln_concentrations)
        opt_dgs = RT * ln_concentrations + self.dG0_f_prime
        return opt_dgs, concentrations, opt_val


class KeggPathway(Pathway):
    
    def __init__(self, S, rids, fluxes, cids, formation_energies=None,
                 reaction_energies=None, cid2bounds=None, c_range=None):
        Pathway.__init__(self, S, formation_energies=formation_energies,
                         reaction_energies=reaction_energies, fluxes=fluxes)
        self.rids = rids
        self.cids = cids
        if cid2bounds:
            self.bounds = [cid2bounds.get(cid, (None, None)) for cid in self.cids]
        else:
            self.bounds = None
        self.cid2bounds = cid2bounds
        self.c_range = c_range
        self.kegg = Kegg.getInstance()

    def GetConcentrationBounds(self, cid):
        lb, ub = None, None
        if cid in self.cid2bounds:
            lb, ub = self.cid2bounds[cid]
        lb = lb or self.c_range[0]
        ub = ub or self.c_range[1]
        return lb, ub

    def GetReactionString(self, r, show_cids=False):
        rid = self.rids[r]
        sparse = dict([(self.cids[c], self.S[r, c]) for c in self.S[r, :].nonzero()[0]])
        if self.fluxes[r] >= 0:
            direction = '=>'
        else:
            direction = '<='
        reaction = Reaction(self.kegg.rid2string(rid), sparse, rid=rid,
                            direction=direction)
        return reaction.to_hypertext(show_cids=show_cids)

    def GetTotalReactionString(self, show_cids=False):
        total_S = np.dot(self.S.T, self.fluxes)
        sparse = dict([(self.cids[c], total_S[c]) for c in total_S.nonzero()[0]])
        reaction = Reaction("Total", sparse, direction="=>")
        return reaction.to_hypertext(show_cids=show_cids)

    def FindMtdf(self, normalization=DeltaGNormalization.DEFAULT):
        """Find the MTDF (Maximal Thermodynamic Driving Force).
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mtdf).
        """
        return self._FindMtdf(self.c_range, self.bounds, normalization)

    def GetTotalReactionEnergy(self):
        return self._GetTotalReactionEnergy(self.c_range, self.bounds)
    
    def GetMaxReactionEnergy(self, min_driving_force=0):
        """
            Maximizes the total pathway dG' (i.e. minimize energetic cost).
            Arguments:
                min_driving_force - the lower limit on each reaction's driving force
                                    (it is common to provide the optimize driving force
                                    in order to find the concentrations that minimize the
                                    cost, without affecting the MTDF).
        """
        return self._GetMaxReactionEnergy(self.c_range, self.bounds, min_driving_force)
        
        
    def FindMinimalFeasibleConcentration(self, cid_to_minimize):
        """
            Compute the smallest ratio between two concentrations which makes the pathway feasible.
            All other compounds except these two are constrained by 'bounds' or unconstrained at all.
        min_conc
            Arguments:
                cid - the CID of the compound whose concentration should be minimized
        
            Returns:
                dGs, concentrations, target-concentration
        """
        index = self.cids.index(cid_to_minimize)
        return Pathway.FindMinimalFeasibleConcentration(self, index,
                                                    self.bounds, self.c_range)

    @staticmethod
    def _EnergyToString(dG):
        if np.isnan(dG):
            return "N/A"
        else:
            return "%.1f" % dG

    @staticmethod
    def _AddProfileToFigure(figure, dGs, fluxes, style, label):
        Nr = dGs.shape[0]
        cum_dG = np.cumsum([0] + [dGs[r, 0] * f for r, f in enumerate(fluxes)])
        plt.plot(np.arange(0.5, Nr + 1), cum_dG, style,
                 figure=figure, label=label)
    
    def PlotProfile(self, concentrations, figure=None):
        if figure is None:
            figure = plt.figure(figsize=(6,6), dpi=100)
        plt.title(r'Thermodynamic profile', figure=figure)
        plt.ylabel(r'cumulative $\Delta G_r$ [kJ/mol]', figure=figure)
        plt.xlabel(r'Reaction', figure=figure)
        
        nonzero_reactions = np.nonzero(self.fluxes)[0]
        nonzero_fluxes = [self.fluxes[i] for i in nonzero_reactions]
        
        plt.xticks(np.arange(0, len(nonzero_fluxes)) + 1,
                   [self.kegg.rid2string(self.rids[i]) for i in nonzero_reactions],
                   fontproperties=FontProperties(size=10), rotation=30)

        if np.isnan(self.dG0_r_prime).any():
            return figure
        
        KeggPathway._AddProfileToFigure(figure,
            self.dG0_r_prime[nonzero_reactions, :], nonzero_fluxes, 'm--', "$\Delta_r G^{'\circ}$")
        
        phys_concentrations = self.GetPhysiologicalConcentrations(self.bounds)
        dGm_r_prime = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        KeggPathway._AddProfileToFigure(figure, 
            dGm_r_prime[nonzero_reactions, :], nonzero_fluxes, 'g--', "$\Delta_r G^{'m}$")

        dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(concentrations)
        KeggPathway._AddProfileToFigure(figure, 
            dG_r_prime[nonzero_reactions, :], nonzero_fluxes, 'b-', "$\Delta_r G^'$")

        plt.legend(loc='lower left')
        return figure
        
    def PlotConcentrations(self, concentrations, figure=None):
        if figure is None:
            figure = plt.figure()
        plt.xscale('log', figure=figure)
        plt.ylabel('Compound KEGG ID', figure=figure)
        plt.xlabel('Concentration [M]', figure=figure)
        plt.yticks(range(self.Nc, 0, -1),
                     ["C%05d" % cid for cid in self.cids],
                     fontproperties=FontProperties(size=8))
        plt.plot(concentrations, range(self.Nc, 0, -1), '*b', figure=figure)

        x_min = concentrations.min() / 10
        x_max = concentrations.max() * 10
        y_min = 0
        y_max = self.Nc + 1
        
        for c, cid in enumerate(self.cids):
            plt.text(concentrations[c, 0] * 1.1, self.Nc - c, self.kegg.cid2name(cid), \
                       figure=figure, fontsize=6, rotation=0)
            b_low, b_up = self.GetConcentrationBounds(cid)
            plt.plot([b_low, b_up], [self.Nc - c, self.Nc - c], '-k', linewidth=0.4)

        if self.c_range is not None:
            plt.axvspan(self.c_range[0], self.c_range[1],
                          facecolor='r', alpha=0.3, figure=figure)
        plt.axis([x_min, x_max, y_min, y_max], figure=figure)
        return figure

    def WriteResultsToHtmlTables(self, html_writer, concentrations):
        self.WriteConcentrationsToHtmlTable(html_writer, concentrations)
        self.WriteProfileToHtmlTables(html_writer, concentrations)

    def WriteConcentrationsToHtmlTable(self, html_writer, concentrations=None):
        #html_writer.write('<b>Compound Concentrations</b></br>\n')
        dict_list = []
        for c, cid in enumerate(self.cids):
            d = {}
            d['KEGG CID'] = '<a href="%s">C%05d</a>' % (self.kegg.cid2link(cid), cid)
            d['Compound Name'] = self.kegg.cid2name(cid)
            lb, ub = self.GetConcentrationBounds(cid)
            d['Concentration LB [M]'] = '%.2e' % lb
            if concentrations is not None:
                d['Concentration [M]'] = '%.2e' % concentrations[c, 0]
            d['Concentration UB [M]'] = '%.2e' % ub
            dict_list.append(d)
        if concentrations is not None:
            headers = ['KEGG CID', 'Compound Name', 'Concentration LB [M]',
                         'Concentration [M]', 'Concentration UB [M]']
        else:
            headers = ['KEGG CID', 'Compound Name', 'Concentration LB [M]',
                       'Concentration UB [M]']
        
        html_writer.write_table(dict_list, headers=headers)
    
    def WriteProfileToHtmlTable(self, html_writer, concentrations=None):
        #html_writer.write('<b>Biochemical Reaction Energies</b></br>\n')
        phys_concentrations = np.ones((len(self.cids), 1)) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if 1 in self.cids:
            # C00001 (water) is an exception, its concentration is always set to 1
            phys_concentrations[self.cids.index(1), 0] = 1.0
        
        dG_r_prime_c = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        headers=["reaction", 'formula', 'flux', 
                 "&#x394;<sub>r</sub>G'<sup>c</sup> [kJ/mol] (%g M)" % self.DEFAULT_PHYSIOLOGICAL_CONC] 
        if concentrations is not None:
            dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(concentrations)
            headers.append("&#x394;<sub>r</sub>G' [kJ/mol]")
        
        dict_list = []
        for r, rid in enumerate(self.rids):
            d = {}
            if type(rid) == types.IntType:
                d['reaction'] = '<a href="%s" title="%s">%s</a>' % \
                            (self.kegg.rid2link(rid),
                             self.kegg.rid2name(rid),
                             self.kegg.rid2string(rid))
            else:
                d['reaction'] = self.kegg.rid2string(rid)
            d['flux'] = "%g" % abs(self.fluxes[r])
            d['formula'] = self.GetReactionString(r, show_cids=False)
            d[headers[3]] = \
                KeggPathway._EnergyToString(np.sign(self.fluxes[r]) * dG_r_prime_c[r, 0])
            if concentrations is not None:
                d[headers[4]] = KeggPathway._EnergyToString(
                        np.sign(self.fluxes[r]) * dG_r_prime[r, 0])
            dict_list.append(d)

        total_dG0 = KeggPathway._EnergyToString(float(np.dot(self.dG0_r_prime.T, self.fluxes)))
        d = {'reaction':'Total',
             'flux':'1',
             'formula':self.GetTotalReactionString(show_cids=False),
             headers[3]:total_dG0}
        if concentrations is not None:
            total_dG = KeggPathway._EnergyToString(float(np.dot(dG_r_prime.T, self.fluxes)))
            d[headers[4]] = total_dG
        dict_list.append(d)
        html_writer.write_table(dict_list, headers=headers)
        
if __name__ == '__main__':
    S = np.array([[-1, 1, 0, 0], [0, -1, 1, 0], [0, 0, 1, -1]])
    dGs = np.array([0, 10, 12, 2]).T
    fluxes = np.array([1, 1, -1])
    rids = [1, 2, 3]
    cids = [1, 2, 3, 4]
    #keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    #dGf, concentrations, protein_cost = keggpath.FindKineticOptimum()
    #print 'protein cost: %g protein units' % protein_cost
    #print 'concentrations: ', concentrations
    #print 'sum(concs): ', sum(concentrations), 'M'
    
    keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    dGf, concentrations, mtdf = keggpath.FindMtdf()
    print 'MDTF: %g' % mtdf

    
    S = np.array([[-1, 1, 0], [0, -1, 1]])
    cvxmod.randseed()
    dGs = np.array(cvxmod.randn(3, 1, std=1000))
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
