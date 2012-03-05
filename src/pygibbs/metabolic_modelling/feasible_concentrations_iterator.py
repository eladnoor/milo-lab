#!/usr/bin/python

import numpy as np

from pygibbs.metabolic_modelling import concentration_optimizer
from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.thermodynamic_constants import default_T, R

RT = R * default_T


class FeasibleConcentrationsIterator(object):
    
    def __init__(self, pathway_model, thermodynamic_data,
                 concentration_bounds):
        self._model = pathway_model
        self._thermo = thermodynamic_data
        self._concentration_bounds = concentration_bounds
        self.S = pathway_model.GetStoichiometricMatrix()
        self.Ncompounds, self.Nrxns = self.S.shape
    
        self.dG0_r_prime = thermodynamic_data.GetDGrTagZero_ForModel(
                self._model)
    
    def __iter__(self):
        mtdf_opt = mtdf_optimizer.MTDFOptimizer(
            self._model, self._thermo)
        res = mtdf_opt.FindMTDF(
            concentration_bounds=self._concentration_bounds)
        
        # Bail entirely if the pathway is infeasible.
        status = res.status
        if status.IsInfeasible():
            return
        # Only return data on successful optimization
        if status.IsSuccessful():
            yield np.matrix(res.ln_concentrations)
        
        conc_opt = concentration_optimizer.ConcentrationOptimizer(
            self._model, self._thermo)
        for i in xrange(self.Ncompounds):
            res = conc_opt.MinimizeConcentration(
                i, concentration_bounds=self._concentration_bounds)
            status = res.status
            if status.IsSuccessful():
                yield np.matrix(res.ln_concentrations)
    
    def Feasible(self, concentrations):
        """Returns True if these concentrations are feasible."""
        dgtag = self.dG0_r_prime + RT * concentrations * self.S
        return (dgtag <= 0).all()