#!/usr/bin/python

import numpy as np

class BaseKineticData(object):
    """Base class for all kinetic data containers."""
    
    def GetKcat(self, reaction_id):
        """Returns the Kcat for this reaction.
        
        Args:
            reaction_id: the identifier of the reaction.
        
        Returns:
            The Kcat (1/s).
        """
        raise NotImplementedError
    
    def GetKm(self, reaction_id, compound_id):
        """Returns the Km for this substrate.
        
        Args:
            reaction_id: the identifier of the reaction.
            compound_id: the identifier of the compound.
            
        Returns:
            The Km (Molar).
        """
        raise NotImplementedError
    
    def GetKcats(self, reaction_ids):
        """Gets the Kcats for all the reaction IDs.
        
        Args:
            reaction_ids: an iterable of reaction identifiers.
        
        Returns:
            A 1xN numpy.matrix of kcats where N = len(reaction_ids).
        """
        N = len(reaction_ids)
        out = np.matrix(np.ones((1, N)))
        for i, rid in enumerate(reaction_ids):
            out[0,i] = self.GetKcat(rid)
        return out
    
    def GetKms(self, reaction_ids, compound_ids):
        """Gets the matrix of Kms for each reaction-compound pair.
        
        Args:
            reaction_ids: an iterable of reaction identifiers.
            compound_ids: an iterable of compound identifiers.
        
        Returns:
            An MxN numpy matrix (M=num compounds, N=num reactions).
        """
        M = len(compound_ids)
        N = len(reaction_ids)
        out = np.matrix(np.ones((M,N)))
        for i, cid in enumerate(compound_ids):
            for j, rid in enumerate(reaction_ids):
                out[i,j] = self.GetKm(rid, cid)
        
        return out
    
    def GetKcatsForModel(self, stoich_model):
        """Returns the Kcat matrix for the given model.
        
        Args:
            stoich_model: a StoichiometricModel object.
        
        Returns:
            A 1xN numpy.matrix of kcats where N = number of reactions.
        """
        reaction_ids = stoich_model.GetReactionIDs()
        return self.GetKcats(reaction_ids)
    
    def GetKmsForModel(self, stoich_model):
        """Returns the Kcat matrix for the given model.
        
        Args:
            stoich_model: a StoichiometricModel object.
        
        Returns:
            A 1xN numpy.matrix of kcats where N = number of reactions.
        """
        reaction_ids = stoich_model.GetReactionIDs()
        compound_ids = stoich_model.GetCompoundIDs()
        return self.GetKms(reaction_ids, compound_ids)
    

class UniformKineticData(BaseKineticData):
    """All enzymes have the same (MM) kinetics."""
    
    def __init__(self, kcat=100, km=1e-4):
        """Initialize the UniformKineticData container.
        
        Args:
            kcat: the global turnover number (1/s).
            km: the global michaelis constant to use (Molar).
        """
        self.kcat = kcat
        self.km = km
    
    def GetKcat(self, reaction_id):
        return self.kcat
    
    def GetKm(self, reaction_id, compound_id):
        return self.km
    
    def GetKcats(self, reaction_ids):
        N = len(reaction_ids)
        return np.matrix(np.ones((1, N))) * self.kcat
    
    def GetKms(self, reaction_ids, compound_ids):
        M = len(compound_ids)
        N = len(reaction_ids)
        return np.matrix(np.ones((M,N))) * self.km
    