#!/usr/bin/python

import logging
import numpy as np

import scipy.optimize as opt

from pygibbs.metabolic_modelling import bounds
from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import general_functors
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.metabolic_modelling import protein_optimized_pathway
from pygibbs.metabolic_modelling import protein_cost_functors


class UnsolvableConvexProblemException(Exception):
    def __init__(self, msg, problem):
        Exception.__init__(self, msg)
        self.problem = problem


class ProteinOptimizer(object):
    """Finds the pathway optimum."""

    DEFAULT_CONC_LB = 1e-6
    DEFAULT_CONC_UB = 1e-2

    def __init__(self, pathway_model, thermodynamic_data, kinetic_data,
                 protein_cost_type=protein_cost_functors.ProteinCostFunc):
        """Initialize the MTDFOptimizer class.
        
        Args:
            pathway_model: the PathwayModel object.
            thermodynamic_data: the ThermodynamicData object.
            kinetic_data: the KineticData object.
            protein_cost_type: the functor type for calculating
              protein cost.
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
        self.protein_cost_type = protein_cost_type
        
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
        
        x0 = initial_concentrations
        if x0 is None:
            x0 = self.GetInitialConditions(my_bounds)
        if x0 is None:
            status = optimized_pathway.OptimizationStatus.Infeasible(
                'MTDF could not be found.')
            return protein_optimized_pathway.ProteinCostOptimizedPathway(
                self._model, self._thermo, my_bounds,
                kinetic_data=self._kinetic_data,
                optimization_status=status)
        
        # Fetch kinetic data
        kcat = self._kinetic_data.GetKcatsForModel(self._model)
        km = self._kinetic_data.GetKmsForModel(self._model)
        masses = self._kinetic_data.GetMassesForModel(self._model)
                
        # Separate out the fixed and non-fixed variables.
        injector = general_functors.FixedVariableInjector(lb, ub, x0)
        initial_conds = injector.GetVariableInitialConds()
        lower_bounds  = injector.GetVariableLowerBounds()
        upper_bounds  = injector.GetVariableUpperBounds()
        
        # Make constraint functions and check feasibility of 
        # initial point.        
        minus_dg_func = general_functors.MinusDG(
            self.S, self.dG0_r_prime, injector, max_dGr=0.0)
        assert (minus_dg_func(initial_conds) >= 0).all()
        
        bounds_func = general_functors.BoundDiffs(
            lower_bounds, upper_bounds)
        assert (bounds_func(initial_conds) >= 0).all()
        
        checkers = [minus_dg_func, bounds_func]
        f_ieq = general_functors.MultiFunctionWrapper(checkers)
        assert (f_ieq(initial_conds) >= 0).all()
        
        # Make the goal function and optimize
        optimization_func = self.protein_cost_type(
            self.S, self.dG0_r_prime, self.fluxes, kcat,
            km, masses, injector)
        initial_conds = np.array(initial_conds.flat)
        initial_func_value = optimization_func(initial_conds)
        initial_deriv = optimization_func.Derivatives(initial_conds)
        approx_deriv = opt.approx_fprime(initial_conds, optimization_func, 1e-8)
        initial_dg = optimization_func.GetDGTag(injector(initial_conds))        
        
        logging.debug('Initial dG %s', initial_dg)
        logging.debug('Initial derivative %s', initial_deriv)
        logging.debug('Approx derivative %s', approx_deriv)
        logging.debug('Diff %s', np.array(initial_deriv) - np.array(approx_deriv))
        logging.debug('Initial optimization value: %.2g', initial_func_value)
        res = opt.fmin_slsqp(optimization_func, initial_conds, 
                             f_ieqcons=f_ieq,
                             full_output=1,
                             fprime=optimization_func.Derivatives,
                             iprint=0,
                             acc=1e-2,
                             iter=1000)
        
        # Failed to optimize.
        status, status_str = res[3:]
        if status != 0:
            status = optimized_pathway.OptimizationStatus.GeneralFailure(status_str)
            return protein_optimized_pathway.ProteinCostOptimizedPathway(
                self._model, self._thermo, my_bounds,
                kinetic_data=self._kinetic_data,
                optimization_status=status)
            
        ln_conc, optimum = res[:2]
        final_func_value = optimization_func(ln_conc)
        final_constraints = (f_ieq(ln_conc) >= 0).all()
        logging.debug('Optimum meets constraints %s', final_constraints)
        logging.debug('Final optimization value: %.2g', final_func_value)
        
        enzyme_levels = optimization_func.GetEnzymeLevels(ln_conc)
        ln_conc = injector(ln_conc)
        stoich_factors = optimization_func.StoichFactor()
        max_rate_factors = optimization_func.MaximalRateFactor()
        kinetic_factors = optimization_func.KineticFactor(ln_conc)
        thermo_factors = optimization_func.ThermoFactor(
            optimization_func.GetDGTag(ln_conc))
        
        return protein_optimized_pathway.ProteinCostOptimizedPathway(
            self._model, self._thermo, my_bounds,
            optimal_value=optimum,
            optimal_ln_metabolite_concentrations=ln_conc,
            protein_levels=enzyme_levels,
            stoich_factors=stoich_factors,
            max_rate_factors=max_rate_factors,
            kinetic_factors=kinetic_factors,
            thermo_factors=thermo_factors,
            kinetic_data=self._kinetic_data)
        
        