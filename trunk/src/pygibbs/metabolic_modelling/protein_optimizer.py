#!/usr/bin/python

import logging
import numpy as np
import pylab

import scipy.optimize as opt

from os import path
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
    
    def __init__(self, *args, **kwargs):
        self.protein_levels = kwargs.pop('protein_levels', None)
        self.max_rate_factors = kwargs.pop('max_rate_factors', None)
        self.kinetic_factors = kwargs.pop('kinetic_factors', None)
        self.thermo_factors = kwargs.pop('thermo_factors', None)        

        optimized_pathway.OptimizedPathway.__init__(self, *args,
                                                    **kwargs)
        
        self.protein_level_filename = '%s_protein_levels.png' % self.slug_name 
    
    def WriteProteinProfile(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        if (self.protein_levels is None or
            self.max_rate_factors is None or 
            self.kinetic_factors is None or
            self.thermo_factors is None):
            self.protein_level_filename = None
            return
        
        pylab.figure()
        
        pylab.yscale('log')
        rxn_range = pylab.arange(len(self.reaction_ids))
        flat_max_rate = np.array(self.max_rate_factors.flat)
        flat_kinetics = np.array(self.kinetic_factors.flat)
        flat_thermo = np.array(self.thermo_factors.flat)
        
        first_height = flat_max_rate
        second_height = flat_max_rate * flat_kinetics
        third_height = flat_max_rate * flat_kinetics * flat_thermo        
        
        bottom_rung = first_height
        second_rung = second_height - first_height
        third_rung = third_height - second_height

        pylab.bar(rxn_range, bottom_rung, color='b',
                  label='Due to maximal rate')
        pylab.bar(rxn_range, second_rung, bottom=first_height, color='g',
                  label='Kinetic penalty')
        pylab.bar(rxn_range, third_rung, bottom=second_height, color='c',
                  label='Thermodynamic penalty')
        
        pylab.xticks(rxn_range + 0.5, self.reaction_ids)
        pylab.xlabel('Reaction Step')
        pylab.ylabel('Protein Units')
        pylab.legend(loc='upper left', prop=LEGEND_FONT)
        
        outfname = path.join(dirname, self.protein_level_filename)
        pylab.savefig(outfname, format='png')
     
    
    def WriteAllGraphs(self, dirname):
        """Writes all graphs to the given directory.
        
        Args:
            dirname: the name of the directory to write to.
        """
        optimized_pathway.OptimizedPathway.WriteAllGraphs(self, dirname)
        self.WriteProteinProfile(dirname)


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
    def Denom(dGr_tag):
        """Calculates the exact denominator."""
        denom = 1 - np.exp(dGr_tag/RT)
        return denom
    
    def MaximalRateFactor(self):
        return self.scaled_fluxes
    
    def ThermoFactor(self, full_x):
        dgtag = self.dG0 + RT * full_x * self.S
        if (dgtag >= 0).any():
            return np.matrix(np.ones(dgtag.shape) * 1e6) 
        return 1.0 / self.Denom(dgtag)
    
    def KineticFactor(self, full_x):
        x_exp = np.matrix(np.exp(full_x))
        scaled_kms = self.km / x_exp.T
        exponentiated = np.power(scaled_kms, self.m_plus)
        prods = np.prod(exponentiated, axis=0)
        return 1 + prods
            
    @staticmethod
    def ApproximateDenom(dGr_tag):
        """Calculates the piecewise linear approximation of the
           denominator of our rate law.
        
        DEPRECATED
        
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
        
        maximal_rate_factor = self.MaximalRateFactor()
        kinetic_factor = self.KineticFactor(my_x)
        thermo_factor = self.ThermoFactor(my_x)
        
        levels = np.multiply(
            np.multiply(maximal_rate_factor, kinetic_factor),
                thermo_factor)
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

    def __init__(self, pathway_model, thermodynamic_data, kinetic_data):
        """Initialize the MTDFOptimizer class.
        
        Args:
            pathway_model: the PathwayModel object.
            thermodynamic_data: the ThermodynamicData object.
            kinetic_data: the KineticData object.
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
        
        # Fetch kinetic data
        kcat = self._kinetic_data.GetKcatsForModel(self._model)
        km = self._kinetic_data.GetKmsForModel(self._model)
                
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
        
        
        enzyme_levels = optimization_func.GetEnzymeLevels(ln_conc)
        ln_conc = injector(ln_conc)
        max_rate_factors = optimization_func.MaximalRateFactor()
        kinetic_factors = optimization_func.KineticFactor(ln_conc)
        thermo_factors = optimization_func.ThermoFactor(ln_conc)
        
        return ProteinCostOptimizedPathway(
            self._model, self._thermo, my_bounds,
            optimal_value=optimum,
            optimal_ln_metabolite_concentrations=ln_conc,
            protein_levels=enzyme_levels,
            max_rate_factors=max_rate_factors,
            kinetic_factors=kinetic_factors,
            thermo_factors=thermo_factors)
        
        