#!/usr/bin/python

import cvxmod
import numpy as np

from cvxmod import atoms
from pygibbs.metabolic_modelling import bounds
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.thermodynamic_constants import default_T, R

RT = R * default_T


class ConcentrationOptimizedPathway(optimized_pathway.OptimizedPathway):
    """Contains the result of a MTDF optimization."""
    OPTIMIZATION_TYPE = "MetaboliteConcentration"
    OPT_UNITS = "M"
    
        

class ConcentrationOptimizer(object):
    """Finds the pathway MTDF."""

    DEFAULT_CONC_LB = 1e-6
    DEFAULT_CONC_UB = 1e-2

    def __init__(self, pathway_model, thermodynamic_data):
        """Initialize the MTDFOptimizer class.
        
        Args:
            pathway_model: the PathwayModel object.
            thermodynamic_data: the ThermodynamicData object.
        """
        self._model = pathway_model
        self._thermo = thermodynamic_data
        self.S = pathway_model.GetStoichiometricMatrix()
        self.Ncompounds, self.Nrxns = self.S.shape
        self.reactions = pathway_model.GetReactionIDs()
        self.compounds = pathway_model.GetCompoundIDs()
        self.fluxes = pathway_model.GetFluxes()
        self.dG0_r_prime = thermodynamic_data.GetDGrTagZero_ForModel(
                self._model)

    def DefaultConcentrationBounds(self):
        """Default Bounds objects."""
        return bounds.Bounds(default_lb=self.DEFAULT_CONC_LB,
                             default_ub=self.DEFAULT_CONC_UB)

    def _LnConcentrationBounds(self, bounds):
        """Make bounds on concentrations from the Bounds object.
        
        Args:
            bounds: a Bounds objects for concentrations.
        
        Returns:
            A 2-tuple (lower bounds, upper bounds) as cvxmod.matrix objects.
        """
        ln_conc_lb = np.matrix(np.zeros((1, self.Ncompounds)))
        ln_conc_ub = np.matrix(np.zeros((1, self.Ncompounds)))
        
        for i, c in enumerate(self.compounds):
            ln_conc_lb[0, i] = bounds.GetLowerBound(c)
            ln_conc_ub[0, i] = bounds.GetUpperBound(c)
        
        ln_conc_lb = np.log(ln_conc_lb)
        ln_conc_ub = np.log(ln_conc_ub)
        
        return cvxmod.matrix(ln_conc_lb), cvxmod.matrix(ln_conc_ub)

    def MinimizeConcentration(self, metabolite_index,
                              concentration_bounds=None):
        """Finds feasible concentrations minimizing the concentration
           of metabolite i.
        
        Args:
            concentration_bounds: the Bounds objects setting concentration bounds.
        """
        problem = cvxmod.problem()
        my_bounds = concentration_bounds or self.DefaultConcentrationBounds()
        
        # Constrain concentrations
        ln_conc = cvxmod.optvar('lnC', rows=1, cols=self.Ncompounds)
        ln_conc_lb, ln_conc_ub = my_bounds.GetLnBounds(self.compounds)
        ln_conc_lb = cvxmod.matrix(ln_conc_lb)
        ln_conc_ub = cvxmod.matrix(ln_conc_ub)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)
        
        # Make the objective
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium.
        S = cvxmod.matrix(self.S)
        
        for i, flux in enumerate(self.fluxes):
            curr_dg0 = self.dG0_r_prime[0, i]
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(curr_dg0):
                continue
            
            curr_dgr = curr_dg0 + RT * ln_conc * S[:, i]
            if flux == 0:
                problem.constr.append(curr_dgr == 0)
            else:                
                problem.constr.append(curr_dgr <= 0)
        
        my_conc = ln_conc[0, metabolite_index]
        problem.objective = cvxmod.minimize(atoms.exp(my_conc))
        status = problem.solve(quiet=True)
        if status != 'optimal':
            status = optimized_pathway.OptimizationStatus.Infeasible(
                'Pathway infeasible given bounds.')
            return ConcentrationOptimizedPathway(
                self._model, self._thermo,
                my_bounds, optimization_status=status)
                    
        optimum = cvxmod.value(problem)
        opt_ln_conc = np.array(cvxmod.value(ln_conc))
        result = ConcentrationOptimizedPathway(
            self._model, self._thermo,
            my_bounds, optimal_value=optimum,
            optimal_ln_metabolite_concentrations=opt_ln_conc)
        return result
        
        
        