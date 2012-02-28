#!/usr/bin/python

import cvxmod
import pylab
import numpy as np

from pygibbs.kegg import Kegg
from pygibbs import kegg_reaction 
from pygibbs.metabolic_modelling import bounds
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.thermodynamic_constants import default_T, R
from toolbox import util
from os import path
from matplotlib.font_manager import FontProperties

RT = R * default_T

LEGEND_FONT = FontProperties(size=8)


class UnsolvableConvexProblemException(Exception):
    def __init__(self, msg, problem):
        Exception.__init__(self, msg)
        self.problem = problem


class MTDFOptimizedPathway(optimized_pathway.OptimizedPathway):
    """Contains the result of a MTDF optimization."""
    OPTIMIZATION_TYPE = "MTDF"
    OPT_UNITS = "kJ / mol"
    
    def GetForwardFraction(self):
        """Computes the thermodynamic efficiency at the MTDF.
        
        Efficiency = (J+ - J-) / (J- + J+)
        According to the formula of Beard and Qian
        """
        return self.CalcForwardFraction(-self.opt_val)
    forward_fraction = property(GetForwardFraction)
    
    def GetNormalization(self):
        return self._normalization
    def SetNormalization(self, norm):
        self._normalization = norm
    normalization = property(GetNormalization, SetNormalization)
    
        

class MTDFOptimizer(object):
    """Finds the pathway MTDF."""

    DEFAULT_CONC_LB = 1e-6
    DEFAULT_CONC_UB = 1e-2

    class DeltaGNormalization(object):
        """Nested class defining dG normalization setting."""
        # Motivation is having the most uniform entropy production
        TIMES_FLUX = 1
         
        # Motivation is requiring more force in reaction that have more flux
        DIVIDE_BY_FLUX = 2
        
        # Motivation is putting limits on the allowed backward/forward fluxes 
        SIGN_FLUX = 3 
        
        # Default setting
        DEFAULT = SIGN_FLUX
        
        @staticmethod
        def NormalizeDGByFlux(dG, flux, normalization):
            if normalization == MTDFOptimizer.DeltaGNormalization.DIVIDE_BY_FLUX:
                return -dG * (1.0 / flux)
            
            if normalization == MTDFOptimizer.DeltaGNormalization.TIMES_FLUX:
                return -dG * flux
            
            if normalization == MTDFOptimizer.DeltaGNormalization.SIGN_FLUX:
                return -dG * np.sign(flux)
            
            raise ValueError('bad value for normalization method: %d',
                             normalization)

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

    def FindMTDF(self, concentration_bounds=None, normalization=None):
        """Finds the MTDF.
        
        Args:
            bounds: the Bounds objects setting concentration bounds.
        """
        problem = cvxmod.problem()
        bounds = concentration_bounds or self.DefaultConcentrationBounds()
        normalization = normalization or self.DeltaGNormalization.DEFAULT
        
        # Constrain concentrations
        ln_conc = cvxmod.optvar('lnC', rows=1, cols=self.Ncompounds)
        # TODO(flamholz): push this method into the bounds object.
        ln_conc_lb, ln_conc_ub = self._LnConcentrationBounds(bounds)
        problem.constr.append(ln_conc >= ln_conc_lb)
        problem.constr.append(ln_conc <= ln_conc_ub)
        
        # Make the objective
        motive_force_lb = cvxmod.optvar('B', 1)
        my_dG0_r_primes = np.matrix(self.dG0_r_prime)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium.
        S = cvxmod.matrix(self.S)
        
        for i, flux in enumerate(self.fluxes):
            
            curr_dg0 = my_dG0_r_primes[0, i]
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(curr_dg0):
                continue
            
            curr_dgr = curr_dg0 + RT * ln_conc * S[:, i]
            if flux == 0:
                problem.constr.append(curr_dgr == 0)
            else:
                motive_force = self.DeltaGNormalization.NormalizeDGByFlux(
                    curr_dgr, flux, normalization)
                
                problem.constr.append(motive_force >= motive_force_lb)
        
        problem.objective = cvxmod.maximize(motive_force_lb)
        status = problem.solve(quiet=True)
        if status != 'optimal':
            raise UnsolvableConvexProblemException(status, problem)
        
        mtdf = cvxmod.value(motive_force_lb)
        opt_ln_conc = np.array(cvxmod.value(ln_conc))
        result = MTDFOptimizedPathway(self._model, self._thermo,
                                      bounds, mtdf, opt_ln_conc)
        result.SetNormalization(normalization)
        return result
        
        
        