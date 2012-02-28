#!/usr/bin/python

import numpy as np

from pygibbs import thermodynamics
from pygibbs import thermodynamic_errors


class BaseThermoData(object):
    """Base class defining the interface."""
    
    def GetDGfTagZero_ForID(self, compound_id):
        """Returns standard formation energy for the given ID.
        
        Args:
            compound_id: the ID of the compound.
        """
        raise NotImplementedError
    
    def GetDGfTagZero_ForIDs(self, compound_ids):
        """Get formation energies in standard conditions.
        
        Args:
            stoich_model: the stoichiometric model.
            
        Returns:
            The formation energies as a 1xlen(compound_ids) row vector.
        """
        formations = np.zeros((1, len(compound_ids)))
        
        for i, id in enumerate(compound_ids):
            formations[0, i] = self.GetDGfTagZero_ForID(id) 
        
        return formations
    
    def GetDGfTagZero_ForModel(self, stoich_model):
        """Get formation energies in standard conditions.
        
        Args:
            stoich_model: the stoichiometric model.
        """
        compound_ids = stoich_model.GetCompoundIDs()
        return self.GetDGfTagZero_ForIDs(compound_ids)
    
    def GetDGrTagZero_ForID(self, reaction_id):
        """Returns standard reaction energy for the given ID.
        
        Args:
            reaction_id: the ID of the reaction.
        """
        raise NotImplementedError

    def GetDGrTagZero_ForIDs(self, reaction_ids):
        """Get reaction energies in standard conditions.
        
        Default implementation.
        
        Args:
            stoich_model: the stoichiometric model.
        """
        rxn_energies = np.zeros((len(reaction_ids),1))
        
        for i, id in enumerate(reaction_ids):
            rxn_energies[i, 0] = self.GetDGrTagZero_ForID(id) 
        
        return rxn_energies
    
    def GetDGrTagZero_ForModel(self, stoich_model):
        """Get reaction energies in standard conditions.
        
        Default implementation.
        
        Args:
            stoich_model: the stoichiometric model.
        """
        reaction_ids = stoich_model.GetReactionIDs()
        return self.GetDGrTagZero_ForIDs(reaction_ids)
    
    
class FormationBasedThermoData(BaseThermoData):
    
    def __init__(self, formation_energies):
        """Initializes formation-energy based thermo data container.
        
        Args:
            formation_energies: a dictionary mapping ids to
                                standard formation energies.
        """                            
        self.formation_energies = formation_energies
    
    def GetDGfTagZero_ForID(self, compound_id):
        """Returns standard formation energy for the given ID.
        
        Args:
            compound_id: the ID of the compound.
        """
        return self.formation_energies.get(compound_id, np.NAN)
    
    def GetDGrTagZero_ForModel(self, stoich_model):
        S = stoich_model.GetStoichiometricMatrix()
        formation_energies = self.GetDGfTagZero_ForModel(stoich_model)
        return thermodynamics.GetReactionEnergiesFromFormationEnergies(
            S, formation_energies)
    


class ReactionBasedThermoData(BaseThermoData):
    
    def __init__(self, reaction_energies):
        """Initializes formation-energy based thermo data container.
        
        Args:
            reaction_energies: a dictionary mapping ids to
                                standard reaction energies.
        """                            
        self.reaction_energies = reaction_energies
    
    def GetDGrTagZero_ForID(self, reaction_id):
        """Returns standard reaction energy for the given ID.
        
        Args:
            reaction_id: the ID of the reaction.
        """
        return self.reaction_energies.get(reaction_id, np.NAN)
    
    
class WrapperThermoData(ReactionBasedThermoData):
    
    def __init__(self, thermo_instance):
        """Initialize the wrapper.
        
        Args:
            thermo_instance: an instance of pygibbs.Thermodynamics. 
        """
        self.thermo_instance = thermo_instance

    def GetDGrTagZero_ForModel(self, stoich_model):
        S = stoich_model.GetStoichiometricMatrix()
        cids = stoich_model.GetCompoundIDs()
        return self.thermo_instance.GetTransfromedReactionEnergies(S, cids)
        
