#!/usr/bin/python

import pylab
import types
from thermodynamic_constants import R, default_T, correction_function, array_transform
from toolbox.util import log_sum_exp
from pygibbs.thermodynamic_constants import dG0_f_Mg

class PseudoisomerMap(object):
    """A map from pseudoisomers to dG values."""
    
    def __init__(self):
        self.dgs = {}
    
    @staticmethod
    def _MakeGroupVectorKey(groupvector):
        return (groupvector.Hydrogens(),
                groupvector.NetCharge(),
                groupvector.Magnesiums())
    
    @staticmethod
    def _MakeKey(nH, z, nMg):
        return (nH, z, nMg)
    
    def AddAll(self, pmap_dict):
        for key, value in pmap_dict.iteritems():
            self.AddKey(key, value)
    
    def Add(self, nH, z, nMg, dG0):
        """For when you don't have a group vector."""
        key = self._MakeKey(nH, z, nMg)
        self.dgs.setdefault(key, []).append(dG0)
    
    def AddKey(self, key, dG0):
        """For when you don't have a group vector."""
        self.dgs.setdefault(key, []).append(dG0)
    
    def AddGroupVector(self, groupvector, dG0):
        """Add a GroupVector with a dG0 value to the map."""
        key = self._MakeGroupVectorKey(groupvector)
        self.dgs.setdefault(key, []).append(dG0)
    
    def Empty(self):
        """Return true if there are no entries in the map."""
        return len(self.dgs) == 0
    
    def Transform(self, pH, pMg, I, T, most_abundant=False):
        """Transform this set of pseudoisomers to a dG value at
           the specified conditions.
        """
        v_dG0 = []
        v_nH  = []
        v_z   = []
        v_mg  = []
        
        for key, dG0_list in self.dgs.iteritems():
            nH, z, nMg = key
            for dG0 in dG0_list:
                v_dG0.append(dG0)
                v_nH.append(nH)
                v_z.append(z)
                v_mg.append(nMg)

        v_dG0 = pylab.array(v_dG0)
        v_nH  = pylab.array(v_nH)
        v_z   = pylab.array(v_z)
        v_mg  = pylab.array(v_mg)

        if most_abundant:
            return min(dG0 + R*T*correction_function(nH, nMg, z, pH, pMg, I, T))
        else:
            return array_transform(v_dG0, v_nH, v_mg, v_z, pH, pMg, I, T)
    
    def TransformMatrix(self, pH, pMg, I, T):
        """Potentially return multiple results..."""
        if (type(pH) != types.ListType and type(I) != types.ListType):
            return self.Transform(pH, pMg, I, T)
        else:
            if (type(pH) != types.ListType):
                pH = [pH]
            if (type(I) != types.ListType):
                I = [I]
    
            dG0_matrix = pylab.zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_matrix[i, j] = self.Transform(pH[i], pMg, I[i], T)
            
            return dG0_matrix

    def ToMatrix(self):
        res = []
        for ((nH, z, nMg), dG0_list) in self.dgs.iteritems():
            for dG0 in dG0_list:
                res.append((nH, z, nMg, dG0))
        return sorted(res)
    
    def __str__(self):
        s = ""
        for ((nH, z, nMg), dG0_list) in self.dgs.iteritems():
            s += "%2d %2d %2d %s\n" % (nH, z, nMg, str(dG0_list))
        return s
    
    def GetdG0(self, nH, z, nMg, T=default_T):
        if (nH, z, nMg) not in self.dgs:
            return None
        
        dG0_list = self.dgs[nH, z, nMg]
        return -(R*T) * log_sum_exp([x/(-R*T) for x in dG0_list])
    
    def GetpKa(self, nH, z, nMg, T=default_T):
        dG0_f_deprotonated = self.GetdG0(nH-1, z-1, nMg, T)
        dG0_f_protonated = self.GetdG0(nH, z, nMg, T)
        if not dG0_f_deprotonated or not dG0_f_protonated:
            return None
        else:
            return (dG0_f_deprotonated - dG0_f_protonated) / (R*T*pylab.log(10))
        
    def GetpK_Mg(self, nH, z, nMg, T=default_T):
        dG0_f_without_Mg = self.GetdG0(nH, z-2, nMg-1, T)
        dG0_f_with_Mg = self.GetdG0(nH, z, nMg, T)
        if not dG0_f_without_Mg or not dG0_f_with_Mg:
            return None
        else:
            return (dG0_f_without_Mg + dG0_f_Mg - dG0_f_with_Mg) / (R*T*pylab.log(10))
        