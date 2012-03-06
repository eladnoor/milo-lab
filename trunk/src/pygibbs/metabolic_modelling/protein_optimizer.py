#!/usr/bin/python

import logging
import numpy as np

import scipy.optimize as opt

from pygibbs.metabolic_modelling import bounds
from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.thermodynamic_constants import default_T, R
from matplotlib.font_manager import FontProperties

RT = R * default_T

LEGEND_FONT = FontProperties(size=8)


class UnsolvableConvexProblemException(Exception):
    def __init__(self, msg, problem):
        Exception.__init__(self, msg)
        self.problem = problem


class ProteinCostOptimizedPathway(optimized_pathway.OptimizedPathway):
    """Contains the result of a protein cost optimization."""
    OPTIMIZATION_TYPE = "Protein Cost"
    OPT_UNITS = "Protein Units / Pathway Flux Units"
    

class FixedVariableInjector(object):
    
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


class EnzymeLevelFunc(object):
    """A callable computing optimal enzyme levels."""
    
    def __init__(self, S, dG0, fluxes, kcat, km,
                 injector):
        self.S = S
        self.m_plus = np.abs(np.clip(S, -1000, 0))
        self.dG0 = dG0
        self.fluxes = fluxes
        self.kcat = kcat
        self.km = km
        self.injector = injector
        self.scaled_fluxes = np.matrix(fluxes / kcat)
        self.Nr, self.Nc = self.S.shape
    
    @staticmethod
    def ExactDenom(dGr_tag):
        """Calculates the exact denominator."""
        denom = 1 - np.exp(dGr_tag/RT)
        return denom
    
    @staticmethod
    def ApproximateDenom(dGr_tag):
        """Calculates the piecewise linear approximation of the
           denominator of our rate law.
           
        Args:
            dGr_tag: a 1xM Numpy matrix of transformed reaction energies.
        
        Returns:
            A 1xM Numpy matrix of calculated denominator values.
        """
        linearized_denom = -dGr_tag / (2*RT)
        over_i = np.where(linearized_denom >= 1.0)
        linearized_denom[over_i] = 1.0
        return linearized_denom
    
    def GetEnzymeLevels(self, x):
        """Get the modeled enzyme levels given concentrations."""
        my_x = self.injector(x)
        x_exp = np.matrix(np.exp(my_x))
        
        dgtag = self.dG0 + RT * my_x * self.S
        if (dgtag >= 0).any():
            return 1e6
                
        #linearized_denom = self.ApproximateDenom(dgtag)
        linearized_denom = self.ExactDenom(dgtag)
        
        scaled_kms = self.km / x_exp.T
        exponentiated = np.power(scaled_kms, self.m_plus)
        prods = np.prod(exponentiated, axis=0)
        numer = 1 + prods
        scaled_numers = np.multiply(self.scaled_fluxes, numer)
        
        inverted_denom = 1 / linearized_denom
        levels = np.multiply(scaled_numers, inverted_denom)
        return levels
    
    def __call__(self, x):
        levels = self.GetEnzymeLevels(x)
        return np.sum(levels)

    
class MinusDG(object):
    """A callable checking in thermodynamic requirements are met."""
    
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
        my_x = self.injector(x)
        
        dgtag = self.dG0r + RT * my_x * self.S
        minus_max = -(dgtag - self.max_dGr)
        return minus_max
    
    
class BoundDiffs(object):
    
    def __init__(self, lb, ub):
        self.lb = lb
        self.ub = ub
    
    def __call__(self, x):
        lb_diff = x - self.lb
        ub_diff = self.ub - x
        return np.matrix(np.hstack([lb_diff, ub_diff]))  
        

class MultiFunctionWrapper(object):
    
    def __init__(self, functions):
        """Initialize the multi-functor.
        
        Args:
            functions: a list of callables.
        """
        self.functions = functions
        
    def __call__(self, x):
        outs = [f(x) for f in self.functions]
        out_mat = np.hstack(outs)
        
        # Required to return a 1-d ndarray.
        return np.array(out_mat.flat)


class ProteinOptimizer(object):
    """Finds the pathway optimum."""

    DEFAULT_CONC_LB = 1e-6
    DEFAULT_CONC_UB = 1e-2

    def __init__(self, pathway_model, thermodynamic_data, kinetic_data=None):
        """Initialize the MTDFOptimizer class.
        
        Args:
            pathway_model: the PathwayModel object.
            thermodynamic_data: the ThermodynamicData object.
        """
        self._model = pathway_model
        self._thermo = thermodynamic_data
        self._kinetic_data = kinetic_data
        self.S = pathway_model.GetStoichiometricMatrix()
        self.Ncompounds, self.Nrxns = self.S.shape
        self.reactions = pathway_model.GetReactionIDs()
        self.compounds = pathway_model.GetCompoundIDs()
        self.fluxes = pathway_model.GetFluxes()
        self.dG0_r_prime = thermodynamic_data.GetDGrTagZero_ForModel(
                self._model)
        
        self.mtdf_opt = mtdf_optimizer.MTDFOptimizer(pathway_model,
                                                     thermodynamic_data)

    def DefaultConcentrationBounds(self):
        """Default Bounds objects."""
        return bounds.Bounds(default_lb=self.DEFAULT_CONC_LB,
                             default_ub=self.DEFAULT_CONC_UB)
    
    def GetInitialConditions(self, concentration_bounds):
        """Solves the MTDF problem to get valid initial conditions.
        
        Returns:
            A Numpy matrix of the initial concentrations to use.
        """
        mtdf_result = self.mtdf_opt.FindMTDF(
            concentration_bounds=concentration_bounds)
        mtdf_value = mtdf_result.opt_val
        if mtdf_value < 0:
            return None

        return np.matrix(mtdf_result.ln_concentrations)
    
    def FindOptimum(self, concentration_bounds=None,
                    initial_concentrations=None):
        """Finds the Optimum.
        
        Args:
            concentration_bounds: the Bounds objects setting concentration bounds.
            initial_conditions: a starting point for the optimization. Must be feasible.
        """
        # Concentration bounds
        my_bounds = concentration_bounds or self.DefaultConcentrationBounds()
        lb, ub = my_bounds.GetLnBounds(self.compounds)
        
        x0 = self.GetInitialConditions(my_bounds)
        if x0 is None:
            status = optimized_pathway.OptimizationStatus.Infeasible(
                'MTDF could not be found.')
            return ProteinCostOptimizedPathway(
                self._model, self._thermo, my_bounds,
                optimization_status=status)
        
        # Kinetic data
        # All Kcat are 100 /s
        kcat = np.ones(self.Nrxns) * 100
        # All Km are 100 uM
        km = np.matrix(np.ones((self.Ncompounds, self.Nrxns))) * 1e-4
                
        # Separate out the fixed and non-fixed variables.
        injector = FixedVariableInjector(lb, ub, x0)
        initial_conds = injector.GetVariableInitialConds()
        lower_bounds  = injector.GetVariableLowerBounds()
        upper_bounds  = injector.GetVariableUpperBounds()
        
        # Make constraint functions and check feasibility of 
        # initial point.        
        minus_dg_func = MinusDG(self.S, self.dG0_r_prime,
                                injector, max_dGr=0.0)
        assert (minus_dg_func(initial_conds) >= 0).all()
        
        bounds_func = BoundDiffs(lower_bounds, upper_bounds)
        assert (bounds_func(initial_conds) >= 0).all()
        
        f_ieq = MultiFunctionWrapper([minus_dg_func, bounds_func])
        assert (f_ieq(initial_conds) >= 0).all()
        
        # Make the goal function and optimize
        optimization_func = EnzymeLevelFunc(self.S, self.dG0_r_prime,
                                            self.fluxes, kcat, km,
                                            injector)
        initial_func_value = optimization_func(initial_conds)
        
        logging.debug('Initial optimization value: %.2g', initial_func_value)
        res = opt.fmin_slsqp(optimization_func, initial_conds, 
                             f_ieqcons=f_ieq,
                             full_output=1,
                             iprint=0)
            
        ln_conc, optimum = res[:2]
        final_func_value = optimization_func(ln_conc)
        final_constraints = (f_ieq(ln_conc) >= 0).all()
        logging.debug('Optimum meets constraints %s', final_constraints)
        logging.debug('Final optimization value: %.2g', final_func_value)
        
        ln_conc = injector(ln_conc)
        return ProteinCostOptimizedPathway(
            self._model, self._thermo, my_bounds,
            optimal_value=optimum, optimal_ln_metabolite_concentrations=ln_conc)
        
        