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
TICK_FONT = FontProperties(size=12)
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
    
    def IsFailure(self):
        return self.status == OptimizationStatus.FAILURE
    failure = property(IsFailure)
    
    def IsInfeasible(self):
        return self.status == OptimizationStatus.INFEASIBLE
    infeasible = property(IsInfeasible)
    
    def IsSuccessful(self):
        return self.status == OptimizationStatus.SUCCESSFUL
    success = property(IsSuccessful)
    
    
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
        self.Ncompounds, self.Nreactions = self.S.shape
        self.status = optimization_status
        self.opt_val = optimal_value
        self.ln_concentrations = optimal_ln_metabolite_concentrations
        
        self.dGr0_tag = np.array(thermodynamic_data.GetDGrTagZero_ForModel(
                self.model))
        self.dGr0_tag_list = list(self.dGr0_tag.flatten())
        self.compound_ids = self.model.GetCompoundIDs()
        self.reaction_ids = self.model.GetReactionIDs()
        self.fluxes = self.model.GetFluxes()
        
        self.slug_name = util.slugify(model.name)
        self.pathway_graph_filename = '%s_graph.svg' % self.slug_name 
        self.thermo_profile_filename = '%s_thermo_profile.png' % self.slug_name
        self.conc_profile_filename = '%s_conc_profile.png' % self.slug_name
        self.kegg = Kegg.getInstance()
        
        self.concentrations = None
        self.dGr_tag = None
        self.dGr_tag_list = None
        self.dGr_bio = None
        self.dGr_bio_list = None
        
        if (self.ln_concentrations is not None and
            self.dGr0_tag is not None):
            self.concentrations = np.exp(self.ln_concentrations)
            conc_correction = RT * self.ln_concentrations * self.S
            self.dGr_tag = np.array(self.dGr0_tag + conc_correction)
            self.dGr_tag_list = list(self.dGr_tag.flatten())
            
            bio_concs = self.bounds.GetBoundsWithDefault(self.compound_ids, default=1e-3)
            bio_correction = RT * np.dot(np.log(bio_concs), self.S)
            self.dGr_bio = np.array(self.dGr0_tag + bio_correction)
            self.dGr_bio_list = list(self.dGr_bio.flatten())
        
    @staticmethod
    def CalcForwardFraction(dg_tag):
        """Computes the thermodynamic efficiency at a particular dG.
        
        Efficiency = (J+ - J-) / (J- + J+)
        According to the formula of Beard and Qian
        """
        return np.tanh(-dg_tag / (2*RT))

    def GetOptimizationType(self):
        return self.OPTIMIZATION_TYPE
    optimization_type = property(GetOptimizationType)
    
    def GetOptimizationUnits(self):
        return self.OPT_UNITS
    optimization_units = property(GetOptimizationUnits)

    def WriteConcProfile(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        pylab.figure()
        pylab.xscale('log')
        pylab.ylabel('Compound KEGG ID')
        pylab.xlabel('Concentration [M]')
        pylab.yticks(range(self.Ncompounds, 0, -1),
                     ["C%05d" % cid for cid in self.compound_ids],
                     fontproperties=TICK_FONT)
        flat_concs = np.array(self.concentrations.flat)
        pylab.plot(flat_concs,
                   range(self.Ncompounds, 0, -1), '*b')

        x_min = self.concentrations.min() / 10
        x_max = self.concentrations.max() * 10
        y_min = 0
        y_max = self.Ncompounds + 1
        
        for c, cid in enumerate(self.compound_ids):
            y_val = self.Ncompounds - c
            pylab.text(self.concentrations[0, c] * 1.1,
                       y_val, self.kegg.cid2name(cid),
                       fontsize=6, rotation=0)
            
            b_low = self.bounds.GetLowerBound(cid)
            b_up = self.bounds.GetUpperBound(cid)
            pylab.plot([b_low, b_up], [y_val, y_val],
                       '-k', linewidth=0.4)

        c_range = self.bounds.GetRange()
        if c_range is not None:
            pylab.axvspan(c_range[0], c_range[1],
                          facecolor='r', alpha=0.3)
        pylab.axis([x_min, x_max, y_min, y_max])
        
        outfname = path.join(dirname, self.conc_profile_filename)
        pylab.savefig(outfname, format='png')

    def WriteThermoProfile(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        pylab.figure()
        dg0_profile = np.cumsum([0] + (self.dGr0_tag * self.fluxes).flatten().tolist())
        dgtag_profile = np.cumsum([0] + (self.dGr_tag * self.fluxes).flatten().tolist())
        dgbio_profile = np.cumsum([0] + (self.dGr_bio * self.fluxes).flatten().tolist())
        rxn_range = pylab.arange(len(self.reaction_ids) + 1)
        pylab.plot(rxn_range, dg0_profile, 'k--',
                   linewidth=3, label='Standard Conditions')
        pylab.plot(rxn_range, dgbio_profile, 'k:',
                   linewidth=3, label='Biological Conditions')
        pylab.plot(rxn_range, dgtag_profile, 'k-',
                   linewidth=3, label='Optimized')
        pylab.xticks(rxn_range[:-1] + 0.5, self.reaction_ids,
                     fontproperties=TICK_FONT)
        pylab.xlabel('Reaction step')
        pylab.ylabel('Cumulative dG (kJ/mol)')
        pylab.legend(loc='upper right', prop=LEGEND_FONT)
        #pylab.ylim((-160, 10))
        pylab.xlim((0, len(self.reaction_ids)-1))
        pylab.grid(b=True)
        
        outfname = path.join(dirname, self.thermo_profile_filename)
        pylab.savefig(outfname.replace('png', 'svg'), format='svg')
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
        if not self.status.success:
            return 
        
        self.WriteThermoProfile(dirname)
        self.WritePathwayGraph(dirname)
        self.WriteConcProfile(dirname)
    
    def MakeReaction(self, rid, vec):
        sparse_reaction = {}
        for j, stoich in enumerate(vec):
            if stoich == 0:
                continue
            sparse_reaction[self.compound_ids[j]] = stoich
        return Reaction(rid, sparse_reaction)
    
    def GetReactionObjects(self):
        """Generator for reactions objects in this pathways model.
        
        Yields:
            kegg_reaction.Reaction objects in order defined.
        """
        for i, rid in enumerate(self.reaction_ids):
            col = np.array(self.S[:,i]).flat
            yield self.MakeReaction(rid, col)
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

    def GetReactionDict(self, i, rid):
        dg0_tag = self.dGr0_tag_list[i]
        dg_tag = self.dGr_tag_list[i]
        d = {'id': rid,
             'dGr0_tag': dg0_tag,
             'dGr_tag': dg_tag,
             'flux': self.fluxes[i],
             'thermo_efficiency': self.CalcForwardFraction(dg_tag)}
        return d
    
    def GetNetDGrTag(self):
        return np.sum(self.fluxes * self.dGr_tag)
    net_dg_tag = property(GetNetDGrTag)

    def ReactionDetails(self):
        """Yields dictionaries describing reactions in the model."""        
        dGr0_tags = self.ReactionStandardEnergyList()
        dGr_tags = self.ReactionTransformedEnergyList()
        for i, rid in enumerate(self.reaction_ids):
            yield self.GetReactionDict(i, rid)
    reaction_details = property(ReactionDetails)
    
    def NetReaction(self):
        """Returns the net reaction as a reaction object."""
        flux_mat = np.matrix(self.fluxes).T
        v_total = np.dot(self.S, flux_mat).flat
        #v_total = np.dot(np.matrix(self.fluxes), self.S).flat
        rxn = self.MakeReaction('Net Reaction', v_total)
        return rxn
    net_reaction = property(NetReaction)
        
        