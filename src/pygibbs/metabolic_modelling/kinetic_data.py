#!/usr/bin/python

import numpy as np
import csv

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


class KineticDataWithDefault(BaseKineticData):
    
    def __init__(self, default_kcat=100, default_km=1e-4):
        """Initialize the UniformKineticData container.
        
        Args:
            default_kcat: the global turnover number (1/s).
            default_km: the global michaelis constant to use (Molar).
        """
        self.default_kcat = default_kcat
        self.default_km = default_km
        self.kcats = {}
        self.kms = {}
    
    def SetKcat(self, reaction_id, kcat):
        self.kcats[reaction_id] = kcat
    
    def SetKm(self, reaction_id, compound_id, km):
        self.kms[(reaction_id, compound_id)] = km
    
    @staticmethod
    def FromFiles(kcat_file, km_file,
                  default_kcat=None, default_km=None):
        ret = KineticDataWithDefault()
        
        f = open(kcat_file)
        r = csv.DictReader(f)
        for row in r:
            rid = row['Short Name']
            val = row['KCAT']
            if val:
                ret.SetKcat(rid, float(val))
        f.close()
        
        f = open(km_file)
        r = csv.DictReader(f)
        r = csv.DictReader(f)
        for row in r:
            rid = row['Short Name']
            cid = row['Substrate CID']
            val = row['KM']
            if val:
                ret.SetKm(rid, cid, float(val))
        f.close()
        
        my_kcat = default_kcat
        my_km = default_km
        if not my_kcat:
            my_kcat = np.mean(ret.kcats.values())
        if not my_km:
            my_km = np.mean(ret.kms.values())
        ret.default_kcat = my_kcat
        ret.default_km = my_km
    
        return ret     
        
    def GetKcat(self, reaction_id):
        if reaction_id in self.kcats:
            return self.kcats[reaction_id]
        return self.default_kcat
    
    def GetKm(self, reaction_id, compound_id):
        key = (reaction_id, compound_id) 
        if key in self.kms:
            return self.kms[key]
        return self.default_km
    
    
        
    