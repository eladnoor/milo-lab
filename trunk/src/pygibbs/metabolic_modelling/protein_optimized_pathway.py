
import pylab
import numpy as np

from os import path
from pygibbs.metabolic_modelling import optimized_pathway
from matplotlib.font_manager import FontProperties

LEGEND_FONT = FontProperties(size=8)


class ProteinCostOptimizedPathway(optimized_pathway.OptimizedPathway):
    """Contains the result of a protein cost optimization."""
    OPTIMIZATION_TYPE = "Protein Cost"
    OPT_UNITS = "Grams Protein / Pathway Flux Units"
    
    def __init__(self, *args, **kwargs):
        self.protein_levels = kwargs.pop('protein_levels', None)
        self.stoich_factors = kwargs.pop('stoich_factors', None)
        self.max_rate_factors = kwargs.pop('max_rate_factors', None)
        self.kinetic_factors = kwargs.pop('kinetic_factors', None)
        self.thermo_factors = kwargs.pop('thermo_factors', None)        
        self.kinetic_data = kwargs.pop('kinetic_data', None)

        optimized_pathway.OptimizedPathway.__init__(self, *args,
                                                    **kwargs)
        
        self.protein_level_filename = '%s_protein_levels.png' % self.slug_name
        self.cumulative_levels_filename = '%s_cumulative_levels.png' % self.slug_name
    
    def GetReactionDict(self, i, rid):
        d = optimized_pathway.OptimizedPathway.GetReactionDict(self, i, rid)
        d['kcat'] = self.kinetic_data.GetKcat(rid)
        d['protein_level'] = self.protein_levels[0, i]
        return d
    
    def ProteinLevelsList(self):
        return list(self.protein_levels.flat)
    protein_levels_list = property(ProteinLevelsList)      
    
    def WriteCumulativeLevels(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        if self.protein_levels is None:
            self.cumulative_levels_filename = None
            return
        
        pylab.figure()
        
        pylab.yscale('log')
        cum_levels = list(np.cumsum(self.protein_levels).flat)
        rxn_range = pylab.arange(len(self.reaction_ids))
        pylab.plot(rxn_range, cum_levels, 'b-')
        pylab.xticks(rxn_range + 0.5, self.reaction_ids)
        pylab.ylabel('Cumulative Grams Protein / Pathway Flux Unit')
        pylab.xlabel('Reaction Step')
                
        outfname = path.join(dirname, self.cumulative_levels_filename)
        pylab.savefig(outfname, format='png')
    
    def WriteProteinProfile(self, dirname):
        """Writes the thermodynamic profile graph.
        
        Writes to "dirname/self.thermo_profile_filename".
        
        Args:
            dirname: the name of the directory to write it to.
        """
        if (self.protein_levels is None or
            self.stoich_factors is None or
            self.max_rate_factors is None or 
            self.kinetic_factors is None or
            self.thermo_factors is None):
            self.protein_level_filename = None
            return
        
        pylab.figure()
        
        pylab.yscale('log')
        rxn_range = pylab.arange(len(self.reaction_ids))
        flat_max_rate = np.array(self.max_rate_factors.flat)
        flat_stoich = np.array(self.stoich_factors.flat)
        flat_kinetics = np.array(self.kinetic_factors.flat)
        flat_thermo = np.array(self.thermo_factors.flat)
        
        first_height  = flat_max_rate
        second_height = first_height * flat_stoich
        third_height  = second_height * flat_kinetics
        fourth_height = third_height * flat_thermo       
        total = np.sum(fourth_height)
        mean = np.mean(fourth_height)
        
        
        bottom_rung = first_height
        second_rung = second_height - first_height
        third_rung  = third_height - second_height
        fourth_rung = fourth_height - third_height

        pylab.bar(rxn_range, bottom_rung,
                  color='w', edgecolor='w',
                  label='Due to maximal rate')
        pylab.bar(rxn_range, second_rung, bottom=first_height,
                  color='#FF5D40', edgecolor='w',
                  label='Due to stoichiometry')
        pylab.bar(rxn_range, third_rung, bottom=second_height,
                  color='#37DD6F', edgecolor='w',
                  label='Due to non-saturation')
        pylab.bar(rxn_range, fourth_rung, bottom=third_height,
                  color='#4186D3', edgecolor='w',
                  label='Due to reverse-reaction')
        
        pylab.xticks(rxn_range + 0.5, self.reaction_ids,
                     fontproperties=optimized_pathway.TICK_FONT)
        pylab.xlabel('Reaction Step')
        pylab.ylabel('Grams Protein / Flux Unit (total %.2g, mean %.2g)' % (total, mean))
        #pylab.ylim((100, 1e4))
        pylab.legend(loc='upper left', prop=LEGEND_FONT)
        
        outfname = path.join(dirname, self.protein_level_filename)
        pylab.savefig(outfname.replace('png', 'svg'), format='svg')
        pylab.savefig(outfname, format='png')
    
    def WriteAllGraphs(self, dirname):
        """Writes all graphs to the given directory.
        
        Args:
            dirname: the name of the directory to write to.
        """
        optimized_pathway.OptimizedPathway.WriteAllGraphs(self, dirname)
        self.WriteProteinProfile(dirname)
        self.WriteCumulativeLevels(dirname)
