#!/usr/bin/python

import numpy as np
import pylab


from matplotlib.font_manager import FontProperties
from os import path
from pygibbs.kegg import Kegg
from pygibbs.kegg_reaction import Reaction 
from pygibbs.thermodynamic_constants import default_RT
from toolbox import util

LEGEND_FONT = FontProperties(size=8)
RT = default_RT

class OptimizationStatus(object):
    
    SUCCESSFUL = 0
    FAILURE    = 1 
    INFEASIBLE = 2
    
    STATUS_STRS = {SUCCESSFUL: 'successful',
                   FAILURE: 'failure',
                   INFEASIBLE: 'infeasible'}
    
    def __init__(self, status,
                 explanation=None):
        self.status = status
        self.explanation = explanation
        
    def __str__(self):
        return '%s (%s)' % (OptimizationStatus.STATUS_STRS.get(self.status),
                            self.explanation)
        
    @staticmethod
    def Successful(explanation=None):
        return OptimizationStatus(OptimizationStatus.SUCCESSFUL,
                                  explanation=explanation)
    
    @staticmethod
    def GeneralFailure(explanation=None):
        return OptimizationStatus(OptimizationStatus.FAILURE,
                                  explanation=explanation)
    
    @staticmethod
    def Infeasible(explanation=None):
        return OptimizationStatus(OptimizationStatus.INFEASIBLE,
                                  explanation=explanation)
    
    

class OptimizedPathway(object):
    """The result of a pathway optimization."""
    # A string denoting the type of optimization.
    OPTIMIZATION_TYPE = None
    OPT_UNITS = None

    def __init__(self,
                 model,
                 thermodynamic_data,
                 metabolite_concentration_bounds,
                 optimization_status=OptimizationStatus.Successful(),
                 optimal_value=None,
                 optimal_ln_metabolite_concentrations=None):
        self.model = model
        self.thermo = thermodynamic_data
        self.bounds = metabolite_concentration_bounds
        self.S = model.GetStoichiometricMatrix()
        self.status = optimization_status
        self.opt_val = optimal_value
        self.ln_concentrations = optimal_ln_metabolite_concentrations
        
        self.dGr0_tag = np.array(thermodynamic_data.GetDGrTagZero_ForModel(
                self.model))
        self.compound_ids = self.model.GetCompoundIDs()
        self.reaction_ids = self.model.GetReactionIDs()
        
        # Don't proceed...
        if (self.ln_concentrations is None or
            self.opt_val is None):
            return
        
        self.concentrations = np.exp(self.ln_concentrations)
        conc_correction = RT * self.ln_concentrations * self.S
        self.dGr_tag = np.array(self.dGr0_tag + conc_correction)
        
        bio_concs = self.bounds.GetBoundsWithDefault(self.compound_ids, default=1e-3)        
        bio_correction = RT * np.dot(np.log(bio_concs), self.S)
        self.dGr_bio = np.array(self.dGr0_tag + bio_correction)
    
        self.kegg = Kegg.getInstance()
        
        slug_name = util.slugify(model.name)
        self.pathway_graph_filename = '%s_graph.svg' % slug_name 
        self.thermo_profile_filename = '%s_mtreturn optimizedf.png' % slug_name

    @staticmethod
    def CalcForwardFraction(dg_tag):
        """Computes the thermodynamic efficiency at a particular dG.
        
        Efficiency = (J+ - J-) / (J- + J+)
        According to the formula of Beard and Qian
        """
        term_1 = 1.0 / (1.0 + np.exp(dg_tag/RT))
        term_2 = 1.0 / (1.0 + np.exp(-dg_tag/RT))
        return term_1 - term_2

    def GetOptimizationType(self):
        return self.OPTIMIZATION_TYPE
    optimization_type = property(GetOptimizationType)
    
    def GetOptimizationUnits(self):
        return self.OPT_UNITS
    optimization_units = property(GetOptimizationUnits)

    def WriteThermoProfile(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        pylab.figure()
        dg0_profile = np.cumsum([0] + self.dGr0_tag.flatten().tolist())
        dgtag_profile = np.cumsum([0] + self.dGr_tag.flatten().tolist())
        dgbio_profile = np.cumsum([0] + self.dGr_bio.flatten().tolist())
        rxn_range = pylab.arange(len(self.reaction_ids) + 1)
        pylab.plot(rxn_range, dg0_profile, 'b--',
                   linewidth=2, label='Standard Conditions')
        pylab.plot(rxn_range, dgbio_profile, 'c--',
                   linewidth=2, label='Biological Conditions')
        pylab.plot(rxn_range, dgtag_profile, 'g-',
                   linewidth=2, label='Optimized')
        pylab.xticks(rxn_range[:-1] + 0.5, self.reaction_ids)
        pylab.xlabel('Reaction step')
        pylab.ylabel('Cumulative dG (kJ/mol)')
        pylab.legend(loc='upper right', prop=LEGEND_FONT)
        
        outfname = path.join(dirname, self.thermo_profile_filename)
        pylab.savefig(outfname, format='png')
    
    def WritePathwayGraph(self, dirname):
        """Writes a pathway DOT graph.
        
        Writes to "dirname/self.pathway_graph_filename".
        
        Args:
            dirname: the name of the directory to write to.
        """
        outfname = path.join(dirname, self.pathway_graph_filename)
        gdot = self.model.GetPathwayGraph()
        gdot.write(outfname, prog='dot', format='svg')
    
    def WriteAllGraphs(self, dirname):
        """Writes all graphs to the given directory.
        
        Args:
            dirname: the name of the directory to write to.
        """
        self.WriteThermoProfile(dirname)
        self.WritePathwayGraph(dirname)
        
    def GetReactionObjects(self):
        """Generator for reactions objects in this pathways model.
        
        Yields:
            kegg_reaction.Reaction objects in order defined.
        """
        for i, rid in enumerate(self.reaction_ids):
            col = np.array(self.S[:,i]).flatten()
            sparse_reaction = {}
            for j, stoich in enumerate(col):
                if stoich == 0:
                    continue
                sparse_reaction[self.compound_ids[j]] = stoich            
            yield Reaction(rid, sparse_reaction)
    reaction_objects = property(GetReactionObjects)

    def GetDGrZeroTag(self):
        """Returns the standard dGr values."""
        return self.dGr0_tag
    
    def GetDGrTag(self):
        """Returns the transformed dGr values at the optimum."""
        return self.dGr_tag
    
    def ConcentrationsList(self):
        """Returns a list of optimum concentrations in order."""
        return list(self.concentrations.flat)

    def CompoundNames(self):
        """Presumes compound IDs are from KEGG."""
        kegg_instance = Kegg.getInstance()
        return map(kegg_instance.cid2name, self.compound_ids)

    def CompoundDetails(self):
        """Yields dictionaries describing compounds in the model."""        
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
        """Returns a list of reaction dG0 values."""
        return list(self.dGr0_tag.flatten())
    
    def ReactionTransformedEnergyList(self):
        """Returns a list of reaction dG' values."""
        return list(self.dGr_tag.flatten())

    def ReactionDetails(self):
        """Yields dictionaries describing reactions in the model."""        
        dGr0_tags = self.ReactionStandardEnergyList()
        dGr_tags = self.ReactionTransformedEnergyList()
        for i, id in enumerate(self.reaction_ids):
            dg0_tag = dGr0_tags[i]
            dg_tag = dGr_tags[i]
            d = {'id': id,
                 'dGr0_tag': dg0_tag,
                 'dGr_tag': dg_tag,
                 'thermo_efficiency': self.CalcForwardFraction(dg_tag)}
            yield d
    reaction_details = property(ReactionDetails)
