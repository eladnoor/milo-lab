#!/usr/bin/python
    

class StoichiometricModel(object):
    """A stoichiometric model (of a pathway or metabolic system).
    
    Contains:
        Contains reaction and compound data.
    """
    
    def __init__(self):
        return

    def GetStoichiometricMatrix(self):
        """Returns the stoichiometric matrix."""
        raise NotImplementedError
    
    def GetReactionIDs(self):
        """Returns the reaction IDs (a list of strings)."""
        raise NotImplementedError
    
    def GetCompoundIDs(self):
        """Returns the compound IDs (a list of strings)."""
        raise NotImplementedError
    
    def GetFluxes(self):
        """Returns the compound IDs."""
        raise NotImplementedError
    
        