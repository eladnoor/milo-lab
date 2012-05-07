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
    
    def GetMassPerActiveSite(self, reaction_id):
        """Returns the mass in kDa of the enzyme per active site.
        
        Args:
            reaction_id: the identifier of the reaction.
        
        Returns:
            The mass per active site (kDa).
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
    
    def GetMassesPerActiveSite(self, reaction_ids):
        """Gets the mass / active site for all the reaction IDs.
        
        Args:
            reaction_ids: an iterable of reaction identifiers.
        
        Returns:
            A 1xN numpy.matrix of masses (kDa) where N = len(reaction_ids).
        """
        N = len(reaction_ids)
        out = np.matrix(np.ones((1, N)))
        for i, rid in enumerate(reaction_ids):
            out[0,i] = self.GetMassPerActiveSite(rid)
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

    def GetMassesForModel(self, stoich_model):
        """Gets the mass / active site for the given model.
        
        Args:
            stoich_model: a StoichiometricModel object.
        
        Returns:
            A 1xN numpy.matrix of masses (kDa) where N = number of reactions.
        """
        reaction_ids = stoich_model.GetReactionIDs()
        return self.GetMassesPerActiveSite(reaction_ids)

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
    
    def __init__(self, kcat=100, km=1e-4, mass=40):
        """Initialize the UniformKineticData container.
        
        Args:
            kcat: the global turnover number (1/s).
            km: the global michaelis constant to use (Molar).
            mass_per_active_site: the enzyme mass per active site.
        """
        self.kcat = kcat
        self.km = km
        self.mass = mass
    
    def GetKcat(self, reaction_id):
        return self.kcat
    
    def GetKm(self, reaction_id, compound_id):
        return self.km
    
    def GetMassPerActiveSite(self, reaction_id):
        return self.mass
    
    def GetKcats(self, reaction_ids):
        N = len(reaction_ids)
        return np.matrix(np.ones((1, N))) * self.kcat
    
    def GetKms(self, reaction_ids, compound_ids):
        M = len(compound_ids)
        N = len(reaction_ids)
        return np.matrix(np.ones((M,N))) * self.km
    
    def GetMassesPerActiveSite(self, reaction_ids):
        N = len(reaction_ids)
        return np.matrix(np.ones((1, N))) * self.mass


class KineticDataWithDefault(BaseKineticData):
    
    def __init__(self, default_kcat=100, default_km=1e-4, default_mass=40):
        """Initialize the UniformKineticData container.
        
        Args:
            default_kcat: the global turnover number (1/s).
            default_km: the global michaelis constant to use (Molar).
            default_mass: the default enzyme mass per active site.
        """
        self.default_kcat = default_kcat
        self.default_km = default_km
        self.default_mass = default_mass
        self.kcats = {}
        self.kms = {}
        self.masses = {}
    
    def SetKcat(self, reaction_id, kcat):
        self.kcats[reaction_id] = kcat
    
    def SetKm(self, reaction_id, compound_id, km):
        self.kms[(reaction_id, compound_id)] = km
    
    def SetMass(self, reaction_id, mass):
        self.masses[reaction_id] = mass
        
    def SetDefaultKcat(self, default_kcat):
        self.default_kcat = default_kcat
    
    def SetDefaultKM(self, default_km):
        self.default_km = default_km
        
    def SetDefaultMass(self, default_mass):
        self.default_mass = default_mass
    
    @staticmethod
    def FromArrenFile(filename):
        ret = KineticDataWithDefault()
        
        f = open(filename)
        r = csv.DictReader(f)
    
        for row in r:
            rid   = row['Short Name'].strip()
            param = row['Parameter'].strip()
            value = row['Value'].strip()
            direc = row['Direction'].strip()
            subst = row['Substrate CID'].strip()
            if not value:
                continue
            value = float(value.replace(',', ''))
            if param == 'MW':
                ret.SetMass(rid, value)
            if param == 'kcat' and direc == 'F':
                ret.SetKcat(rid, value)
            if param == 'KM' and subst:
                # Convert KMs to M
                ret.SetKm(rid, subst, value/1000.0)

        my_kcat = np.exp(np.mean(np.log(ret.kcats.values())))
        my_km = np.exp(np.mean(np.log(ret.kms.values())))
        masses = np.array(ret.masses.values())
        masses = masses[np.nonzero(masses)]
        my_mass = np.exp(np.mean(np.log(masses)))
        ret.default_kcat = my_kcat
        ret.default_km   = my_km
        ret.default_mass = my_mass

        f.close()
        return ret
        
    @staticmethod
    def FromFiles(kcat_file, km_file,
                  default_kcat=None,
                  default_km=None,
                  default_mass=None):
        ret = KineticDataWithDefault()
        
        f = open(kcat_file)
        r = csv.DictReader(f)
        for row in r:
            rid = row['Short Name']
            kcat = row['KCAT']
            mass = row['Mass']
            if kcat:
                ret.SetKcat(rid, float(kcat))
            if mass:
                ret.SetMass(rid, float(mass))
        f.close()
        
        f = open(km_file)
        r = csv.DictReader(f)
        for row in r:
            rid = row['Short Name']
            cid = row['Substrate CID']
            val = row['KM']
            if val:
                ret.SetKm(rid, cid, float(val))
        f.close()
        
        my_kcat = default_kcat
        my_km   = default_km
        my_mass = default_mass
        if not my_kcat:
            my_kcat = np.exp(np.mean(np.log(ret.kcats.values())))
        if not my_km:
            my_km = np.exp(np.mean(np.log(ret.kms.values())))
        if not my_mass:
            my_mass = np.exp(np.mean(np.log(ret.masses.values())))
        ret.default_kcat = my_kcat
        ret.default_km   = my_km
        ret.default_mass = my_mass
    
        return ret     
        
    def GetKcat(self, reaction_id):
        return self.kcats.get(reaction_id, self.default_kcat)
    
    def GetKm(self, reaction_id, compound_id):
        key = (reaction_id, compound_id)
        return self.kms.get(key, self.default_km)
    
    def GetMassPerActiveSite(self, reaction_id):
        return self.masses.get(reaction_id, self.default_mass)
    
    
        
    