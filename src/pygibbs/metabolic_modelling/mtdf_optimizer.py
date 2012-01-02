#!/usr/bin/python

import cvxmod
import numpy as np

from pygibbs.kegg import Kegg
from pygibbs.metabolic_modelling import bounds
from pygibbs.thermodynamic_constants import default_T, R

RT = R * default_T


class UnsolvableConvexProblemException(Exception):
    def __init__(self, msg, problem):
        Exception.__init__(self, msg)
        self.problem = problem


class MTDFResult(object):
    """A class containing the result of an MTDF optimization."""
    
    def __init__(self, model, thermodynamic_data,
                 bounds, normalization,
                 mtdf_value, ln_concentrations):
        self.model = model
        self.thermo = thermodynamic_data
        self.bounds = bounds
        self.normalization = normalization
        self.S = model.GetStoichiometricMatrix()
        self.mtdf_value = mtdf_value
        self.ln_concentrations = ln_concentrations
        self.dGr0_tag = thermodynamic_data.GetDGrTagZero_ForModel(
                self.model)
        
        self.concentrations = np.exp(self.ln_concentrations)
        conc_correction = RT * np.dot(self.S, self.ln_concentrations)
        self.dGr_tag = self.dGr0_tag + conc_correction
        
        self.compound_ids = self.model.GetCompoundIDs()
        self.reaction_ids = self.model.GetReactionIDs()
    
    def GetConcentrations(self):
        return self.concentrations
    
    def GetLnConcentrations(self):
        return self.ln_concentrations
    
    def GetMTDF(self):
        return self.mtdf_value
    mtdf = property(GetMTDF)
    
    def GetThermodynamicEfficiency(self):
        """Computes the thermodynamic efficiency at the MTDF.
        
        Efficiency = (J+ - J-) / (J- + J+)
        According to the formula of Beard and Qian
        """
        term_1 = 1.0 / (1.0 + np.exp(-self.mtdf/RT))
        term_2 = 1.0 / (1.0 + np.exp(self.mtdf/RT))
        return term_1 - term_2
    thermodynamic_efficiency = property(GetThermodynamicEfficiency)
    
    def GetDGrZeroTag(self):
        """Returns the standard dGr values."""
        return self.dGr0_tag
    
    def GetDGrTag(self):
        """Returns the transformed dGr values at the optimum."""
        return self.dGr_tag
    
    def ConcentrationsList(self):
        return list(self.concentrations.flatten())
    
    def CompoundNames(self):
        """Presumes compound IDs are from KEGG."""
        kegg_instance = Kegg.getInstance()
        return map(kegg_instance.cid2name, self.compound_ids)
    
    def CompoundDetails(self):
        names = self.CompoundNames()
        concentrations = self.ConcentrationsList()
        for i, id in enumerate(self.compound_ids):
            ub = self.bounds.GetUpperBound(id)
            lb = self.bounds.GetLowerBound(id)
            conc = concentrations[i]
            conc_class = None
            if ub == lb:
                conc_class = 'fixedConc'
            elif abs(ub-conc) < 1e-9:
                conc_class = 'concAtUB'
            elif abs(conc-lb) < 1e-9:
                conc_class = 'concAtLB'
            d = {'id': id,
                 'name': names[i],
                 'concentration': conc,
                 'class': conc_class,
                 'ub': ub,
                 'lb': lb}
            yield d
    compound_details = property(CompoundDetails)

    def ReactionStandardEnergyList(self):
        v = list(self.dGr0_tag.flatten())
        return v
    
    def ReactionTransformedEnergyList(self):
        v = list(self.dGr_tag.flatten())
        return v

    def ReactionDetails(self):
        dGr0_tags = self.ReactionStandardEnergyList()
        dGr_tags = self.ReactionTransformedEnergyList()
        for i, id in enumerate(self.reaction_ids):
            diff = abs(dGr_tags[i] + self.mtdf)
            d = {'id': id,
                 'dGr0_tag': dGr0_tags[i],
                 'dGr_tag': dGr_tags[i],
                 'at_mtdf': diff < 1e-6}
            yield d
    reaction_details = property(ReactionDetails)
        

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
        self.Nr, self.Nc = self.S.shape
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
        ln_conc_lb = []
        ln_conc_ub = []
        
        for c in self.compounds:
            ln_conc_lb.append(bounds.GetLowerBound(c))
            ln_conc_ub.append(bounds.GetUpperBound(c))
        
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
        ln_conc = cvxmod.optvar('lnC', self.Nc)
        # TODO(flamholz): push this method into the bounds object.
        ln_conc_lb, ln_conc_ub = bounds.GetLnBounds(self.compounds)
        ln_conc_lb, ln_conc_ub = cvxmod.matrix(ln_conc_lb), cvxmod.matrix(ln_conc_ub) 
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
                motive_force = self.DeltaGNormalization.NormalizeDGByFlux(
                    curr_dgr, flux, normalization)
                
                problem.constr.append(motive_force >= motive_force_lb)
        
        problem.objective = cvxmod.maximize(motive_force_lb)
        status = problem.solve(quiet=True)
        if status != 'optimal':
            raise UnsolvableConvexProblemException(status, problem)
        
        mtdf = cvxmod.value(motive_force_lb)
        opt_ln_conc = np.array(cvxmod.value(ln_conc))
        return MTDFResult(self._model, self._thermo,
                          bounds, normalization,
                          mtdf, opt_ln_conc)
        
        