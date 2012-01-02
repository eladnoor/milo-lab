#!/usr/bin/python
    
import numpy as np


class StoichiometricModel(object):
    """A stoichiometric model (of a pathway or metabolic system).
    
    Contains: reaction, compound and flux data along with the
              stoichiometric matrix.
    """
    
    def __init__(self, S, reaction_ids, compound_ids,
                 fluxes=None, name=None):
        """Initialize the stoichiometric model.
        
        Args:
            S: the stoichiometrix matrix.
               Reactions are on the rows, compounds on the columns.
            reaction_ids: the ids/names of the reactions (rows).
            compound_ids: the ids/names of the compounds (columns).
            fluxes: the list of relative fluxes through all reactions.
                    if not supplied, assumed to be 1.0 for all reactions.
            name: a string name for this model.
        """
        self.S = S
        self.reaction_ids = reaction_ids
        self.compound_ids = compound_ids
        self.fluxes = fluxes
        self.Nr = len(self.reaction_ids)
        self.Nc = len(self.compound_ids)
        self.name = name
        
        expected_Nr, expected_Nc = self.S.shape
        if self.Nr != expected_Nr:
            raise ValueError('Number of rows does not match number of reactions')
        if self.Nc != expected_Nc:
            raise ValueError('Number of columns does not match number of compounds')
        
        if self.fluxes is None:
            self.fluxes = np.ones((self.Nr, 1)) 

    def GetStoichiometricMatrix(self):
        """Returns the stoichiometric matrix."""
        return self.S
    
    def GetReactionIDs(self):
        """Returns the reaction IDs (a list of strings)."""
        return self.reaction_ids
    
    def GetCompoundIDs(self):
        """Returns the compound IDs (a list of strings)."""
        return self.compound_ids
    
    def GetFluxes(self):
        """Returns the compound IDs."""
        return self.fluxes
    
        