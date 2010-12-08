#!/usr/bin/python

import pylab
import types
import thermodynamics

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
    def _MakeKey(nH, z, mgs):
        return (nH, z, mgs)
    
    def AddAll(self, pmap_dict):
        for key, value in pmap_dict.iteritems():
            self.AddKey(key, value)
    
    def Add(self, nH, z, mgs, dG0):
        """For when you don't have a group vector."""
        key = self._MakeKey(nH, z, mgs)
        self.dgs.setdefault(key, []).append(dG0)
    
    def AddKey(self, key, dG0):
        """For when you don't have a group vector."""
        self.dgs.setdefault(key, []).append(dG0)
    
    def AddGroupVector(self, groupvector, dG0):
        key = self._MakeGroupVectorKey(groupvector)
        self.dgs.setdefault(key, []).append(dG0)
    
    def Empty(self):
        return len(self.dgs) == 0
    
    def Transform(self, pH, pMg, I, T, most_abundant=False):
        v_dG0 = []
        v_nH  = []
        v_z   = []
        v_mg  = []
        
        for key, dG0_list in self.dgs.iteritems():
            nH, z, mgs = key
            for dG0 in dG0_list:
                v_dG0.append(dG0)
                v_nH.append(nH)
                v_z.append(z)
                v_mg.append(mgs)

        v_dG0 = pylab.array(v_dG0)
        v_nH  = pylab.array(v_nH)
        v_z   = pylab.array(v_z)
        v_mg  = pylab.array(v_mg)

        if most_abundant:
            return min(v_dG0 / (-thermodynamics.R*T) +
                       thermodynamics.Thermodynamics.correction_function(v_nH, v_mg, v_z,
                                                                         pH, pMg, I, T))
        else:
            return thermodynamics.Thermodynamics.array_transform(v_dG0, v_nH, v_mg, v_z,
                                                                 pH, pMg, I, T)
    
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
        for ((nH, z, mgs), dG0_list) in self.dgs.iteritems():
            for dG0 in dG0_list:
                res.append((nH, z, mgs, dG0))
        return sorted(res)
    
    def __str__(self):
        s = ""
        for ((nH, z, mgs), dG0_list) in self.dgs.iteritems():
            s += "%2d %2d %2d %s\n" % (nH, z, mgs, str(dG0_list))
        return s