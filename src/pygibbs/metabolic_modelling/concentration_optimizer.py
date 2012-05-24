#!/usr/bin/python

import cvxpy
import numpy as np

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

    def MinimizeConcentration(self, metabolite_index=None,
                              concentration_bounds=None):
        """Finds feasible concentrations minimizing the concentration
           of metabolite i.
        
        Args:
            metabolite_index: the index of the metabolite to minimize.
                if == None, minimize the sum of all concentrations.
            concentration_bounds: the Bounds objects setting concentration bounds.
        """
        my_bounds = concentration_bounds or self.DefaultConcentrationBounds()

        ln_conc = cvxpy.variable(m=1, n=self.Ncompounds, name='lnC')
        ln_conc_lb, ln_conc_ub = my_bounds.GetLnBounds(self.compounds)
        constr = [cvxpy.geq(ln_conc, cvxpy.matrix(ln_conc_lb)),
                  cvxpy.leq(ln_conc, cvxpy.matrix(ln_conc_ub))]
        
        my_dG0_r_primes = np.matrix(self.dG0_r_prime)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium.
        S = np.matrix(self.S)
        
        for i, flux in enumerate(self.fluxes):
            
            curr_dg0 = my_dG0_r_primes[0, i]
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(curr_dg0):
                continue
            
            rcol = cvxpy.matrix(S[:, i])
            curr_dgr = curr_dg0 + RT * ln_conc * rcol
            if flux == 0:
                constr.append(cvxpy.eq(curr_dgr, 0.0))
            else:
                constr.append(cvxpy.leq(curr_dgr, 0.0))        
        
        objective = None
        if metabolite_index is not None:
            my_conc = ln_conc[0, metabolite_index]
            objective = cvxpy.minimize(cvxpy.exp(my_conc))
        else:
            objective = cvxpy.minimize(
                cvxpy.sum(cvxpy.exp(ln_conc)))
        
        name = 'CONC_OPT'
        if metabolite_index:
            name = 'CONC_%d_OPT' % metabolite_index
        problem = cvxpy.program(objective, constr, name=name)
        optimum = problem.solve(quiet=True)        

        """
        status = problem.solve(quiet=True)
        if status != 'optimal':
            status = optimized_pathway.OptimizationStatus.Infeasible(
                'Pathway infeasible given bounds.')
            return ConcentrationOptimizedPathway(
                self._model, self._thermo,
                my_bounds, optimization_status=status)
        """
        opt_ln_conc = np.matrix(np.array(ln_conc.value))
        result = ConcentrationOptimizedPathway(
            self._model, self._thermo,
            my_bounds, optimal_value=optimum,
            optimal_ln_metabolite_concentrations=opt_ln_conc)
        return result
        
        
        