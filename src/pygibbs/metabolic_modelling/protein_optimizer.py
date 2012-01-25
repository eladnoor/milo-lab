#!/usr/bin/python

import pylab
import numpy as np

import scipy.optimize as opt

from pygibbs.kegg import Kegg
from pygibbs import kegg_reaction 
from pygibbs.metabolic_modelling import bounds
from pygibbs.metabolic_modelling import mtdf_optimizer
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


class Result(object):
    """A class containing the result of an MTDF optimization."""
    
    def __init__(self, model, thermodynamic_data,
                 bounds, mtdf_value, ln_concentrations):
        self.model = model
        self.thermo = thermodynamic_data
        self.bounds = bounds
        self.S = model.GetStoichiometricMatrix()
        self.mtdf_value = mtdf_value
        self.ln_concentrations = ln_concentrations
        self.dGr0_tag = thermodynamic_data.GetDGrTagZero_ForModel(
                self.model)

        self.compound_ids = self.model.GetCompoundIDs()
        self.reaction_ids = self.model.GetReactionIDs()
        
        self.concentrations = np.exp(self.ln_concentrations)
        conc_correction = RT * np.dot(self.S, self.ln_concentrations)
        self.dGr_tag = self.dGr0_tag + conc_correction
        
        bio_concs = self.bounds.GetBoundsWithDefault(self.compound_ids, default=1e-3)        
        bio_correction = RT * np.dot(self.S, np.log(bio_concs))
        self.dGr_bio = self.dGr0_tag + bio_correction

        self.kegg = Kegg.getInstance()
        
        slug_name = util.slugify(model.name)
        self.graph_filename = '%s_graph.svg' % slug_name 
        self.mtdf_filename = '%s_mtdf.png' % slug_name
        
    def WriteMTDFGraph(self, dirname):
        pylab.figure()
        dg0_profile = np.cumsum([0] + list(self.dGr0_tag))
        dgtag_profile = np.cumsum([0] + list(self.dGr_tag))
        dgbio_profile = np.cumsum([0] + list(self.dGr_bio))
        rxn_range = pylab.arange(len(self.reaction_ids) + 1)
        pylab.plot(rxn_range, dg0_profile, 'b--',
                   linewidth=2, label='Standard Conditions')
        pylab.plot(rxn_range, dgbio_profile, 'c--',
                   linewidth=2, label='Biological Conditions')
        pylab.plot(rxn_range, dgtag_profile, 'g-',
                   linewidth=2, label='MTDF-optimized')
        pylab.xticks(rxn_range[:-1] + 0.5, self.reaction_ids)
        pylab.xlabel('Reaction step')
        pylab.ylabel('Cumulative dG (kJ/mol)')
        pylab.legend(loc='upper right', prop=LEGEND_FONT)
        
        outfname = path.join(dirname, self.mtdf_filename)
        pylab.savefig(outfname, format='png')
        
    def WritePathwayGraph(self, dirname):
        outfname = path.join(dirname, self.graph_filename)
        gdot = self.model.GetPathwayGraph()
        gdot.write(outfname, prog='dot', format='svg')
        
    def WriteAllGraphs(self, dirname):
        self.WriteMTDFGraph(dirname)
        self.WritePathwayGraph(dirname)
    
    def GetReactionObjects(self):
        for i, rid in enumerate(self.reaction_ids):
            row = self.S[i,:].flatten()
            sparse_reaction = {}
            for j, stoich in enumerate(row):
                if stoich == 0:
                    continue
                sparse_reaction[self.compound_ids[j]] = stoich            
            yield kegg_reaction.Reaction(rid, sparse_reaction)
    reaction_objects = property(GetReactionObjects)
    
    def GetConcentrations(self):
        return self.concentrations
    
    def GetLnConcentrations(self):
        return self.ln_concentrations
    
    def GetMTDF(self):
        return self.mtdf_value
    mtdf = property(GetMTDF)
    
    def GetForwardFraction(self):
        """Computes the thermodynamic efficiency at the MTDF.
        
        Efficiency = (J+ - J-) / (J- + J+)
        According to the formula of Beard and Qian
        """
        term_1 = 1.0 / (1.0 + np.exp(-self.mtdf/RT))
        term_2 = 1.0 / (1.0 + np.exp(self.mtdf/RT))
        return term_1 - term_2
    forward_fraction = property(GetForwardFraction)
    
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


class FixedVariableInjector(object):
    
    def __init__(self, variable_i, fixed_i, fixed_vals):
        self.variable_i = np.array(variable_i)
        self.fixed_i = np.array(fixed_i)
        self.fixed_vals = np.array(fixed_vals)
        self.n = len(self.fixed_i) + len(self.variable_i)
    
    def __call__(self, x):
        out = np.zeros(self.n)
        if self.fixed_i.any():
            out[self.fixed_i] = self.fixed_vals
        out[self.variable_i] = x
        return out
    

class EnzymeLevelFunc(object):
    """A callable computing optimal enzyme levels."""
    
    def __init__(self, S, dG0, fluxes, kcat, km,
                 fixed_var_injector):
        self.S = S
        self.dG0 = dG0
        self.fluxes = fluxes
        self.kcat = kcat
        self.km = km
        self.injector = fixed_var_injector
        self.scaled_fluxes = fluxes / kcat
        self.Nr, self.Nc = self.S.shape
    
    def __call__(self, x):
        my_x = self.injector(x)
        x_col = my_x.reshape(self.Nc, 1)
        x_exp = np.exp(my_x)
        
        dgtag = self.dG0 + RT * np.dot(self.S, x_col)
        if (dgtag >= 0).any():
            return 1e6
        
        linearized_denom = -dgtag / RT
        over_i = np.where(linearized_denom > 1.0)[0]
        linearized_denom[over_i] = 1.0
        
        numer = np.ones(self.Nr)
        for i in xrange(self.Nr):
            for j, stoich in enumerate(self.S[i, :]):
                if stoich >= 0:
                    continue
                
                abs_stoich = np.abs(stoich)
                scaled_km = self.km[i,j] / x_exp[j]
                numer[i] = numer[i] * (scaled_km**abs_stoich)
        
        numer = 1 + numer
        enzyme_levels = self.scaled_fluxes.copy()
        for i in xrange(len(self.dG0)):
            level_i = enzyme_levels[i] * (numer[i] / linearized_denom[i])
            enzyme_levels[i] = level_i
            
        return np.sum(enzyme_levels)
    
    
class MinusDG(object):
    """A callable checking in thermodynamic requirements are met."""
    
    def __init__(self, S, dG0,
                 fixed_var_injector,
                 max_dG=0.0):
        self.S = S
        self.injector = fixed_var_injector
        self.Nr, self.Nc = self.S.shape
        self.dG0 = dG0
        self.max_dG = max_dG
    
    def __call__(self, x):
        my_x = self.injector(x)
        x_col = my_x.reshape((self.Nc, 1))
        dgtag = self.dG0 + RT * np.dot(self.S, x_col)
        minus_max = -(dgtag - self.max_dG)
        return minus_max.flatten()
    
    
class BoundDiffs(object):
    
    def __init__(self, lb, ub):
        self.lb = np.array(lb)
        self.ub = np.array(ub)
    
    def __call__(self, x):
        lb_diff = x - self.lb
        ub_diff = self.ub - x
        return np.hstack([lb_diff, ub_diff])  
        

class MultiFunctionWrapper(object):
    
    def __init__(self, functions):
        self.functions = functions
        
    def __call__(self, x):
        outs = [f(x) for f in self.functions]
        return np.hstack(outs)


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
        self.Nr, self.Nc = self.S.shape
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

    def _GetLnConcentrationBounds(self, bounds):
        """Make bounds on concentrations from the Bounds object.
        
        Args:
            bounds: a Bounds objects for concentrations.
        
        Returns:
            A 2-tuple (lower bounds, upper bounds) as objects.
        """
        conc_lb = []
        conc_ub = []
        
        for c in self.compounds:
            conc_lb.append(bounds.GetLowerBound(c))
            conc_ub.append(bounds.GetUpperBound(c))
        
        return np.log(conc_lb), np.log(conc_ub)

    def FindOptimum(self, concentration_bounds=None):
        """Finds the Optimum.
        
        Args:
            concentration_bounds: the Bounds objects setting concentration bounds.
        """
        bounds = concentration_bounds or self.DefaultConcentrationBounds()
        
        lb, ub = self._GetLnConcentrationBounds(bounds)
        kcat = np.ones((self.Nr, 1)) * 100        # All Kcat are 100 /s
        km = np.ones((self.Nr, self.Nc)) * 1e-4   # All Km are 100 uM
                
        # Initial solution is are MTDF concentrations
        mtdf_result = self.mtdf_opt.FindMTDF(concentration_bounds)
        x0 = mtdf_result.ln_concentrations.flatten()
        
        # Separate out the fixed and non-fixed variables.
        fixed_i = np.where(ub == lb)[0]
        variable_i = np.where(ub != lb)[0]
        fixed_x = ub[fixed_i]
        initial_ln_concs = x0[variable_i]
        lower_bounds = lb[variable_i]
        upper_bounds = ub[variable_i]
        
        # Fix numeric issues near the bounds.
        lb_problem_i = np.where(initial_ln_concs < lower_bounds)[0]
        ub_problem_i = np.where(initial_ln_concs > upper_bounds)[0]
        initial_ln_concs[lb_problem_i] = lower_bounds[lb_problem_i]
        initial_ln_concs[ub_problem_i] = upper_bounds[ub_problem_i]
        
        # Class that injects fixed variable values
        injector = FixedVariableInjector(variable_i, fixed_i, fixed_x)
        initial_conds = np.array(initial_ln_concs)
        
        # Make constraint functions and check feasibility of 
        # initial point.        
        minus_dg_func = MinusDG(self.S, self.dG0_r_prime,
                                injector, max_dG=0.0)
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
        print 'Initial optimization value: %.2g' % initial_func_value
        res = opt.fmin_slsqp(optimization_func, initial_conds, 
                             f_ieqcons=f_ieq,
                             full_output=1,
                             iprint=2)
        
        ln_conc, optimum = res[:2]
        final_func_value = optimization_func(ln_conc)
        final_constraints = (f_ieq(ln_conc) >= 0).all()
        print 'Optimum meets constraints', final_constraints
        print 'Final optimization value: %.2g' % final_func_value
        
        ln_conc = injector(ln_conc)
        ln_conc = ln_conc.reshape((self.Nc, 1))
        return Result(self._model, self._thermo,
                      bounds, optimum, ln_conc)
        
        