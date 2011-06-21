#!/usr/bin/python

import pylab
import types
from thermodynamic_constants import R, default_T, correction_function
from toolbox.util import log_sum_exp
from pygibbs.thermodynamic_constants import dG0_f_Mg

class PseudoisomerMap(object):
    """A map from pseudoisomers to dG values."""
    
    def __init__(self, nH=None, z=None, nMg=None, dG0=None, ref=''):
        self.dgs = {}
        self.refs = {}
        
        if (nH != None and
            z != None and
            nMg != None and
            dG0 != None):
            self.Add(nH, z, nMg, dG0, ref)
    
    @staticmethod
    def _MakeGroupVectorKey(groupvector):
        return (groupvector.Hydrogens(),
                groupvector.NetCharge(),
                groupvector.Magnesiums())
    
    def GetAllPseudoisomers(self):
        """Returns tuples of (key, dG0) where the key is
           (nH, charge, nMg)
        """
        for k, v in self.dgs.iteritems():
            nH, z, nMg = k
            yield nH, z, nMg, v
    
    all_pseudoisomers = property(GetAllPseudoisomers)
    
    @staticmethod
    def _MakeKey(nH, z, nMg):
        return (nH, z, nMg)
    
    def AddKey(self, key, dG0, ref=None):
        """For when you don't have a group vector."""
        self.dgs.setdefault(key, []).append(dG0)
        self.refs[key] = ref # stores only the ref for the last added dG0
    
    def Add(self, nH, z, nMg, dG0, ref=''):
        """For when you don't have a group vector."""
        key = self._MakeKey(nH, z, nMg)
        self.AddKey(key, dG0, ref)
    
    def AddGroupVector(self, groupvector, dG0, ref=''):
        """Add a GroupVector with a dG0 value to the map."""
        key = self._MakeGroupVectorKey(groupvector)
        self.AddKey(key, dG0, ref)
    
    def Squeeze(self, T=default_T):
        """Groups together all pseudoisomers that have the same key"""
        squeezed_dgs = {}
        for k, dG0_list in self.dgs.iteritems():
            squeezed_dG0 = -(R*T) * log_sum_exp([x/(-R*T) for x in dG0_list])
            squeezed_dgs[k] = [squeezed_dG0] 
        self.dgs = squeezed_dgs
    
    def Empty(self):
        """Return true if there are no entries in the map."""
        return len(self.dgs) == 0
    
    def _Transform(self, pH, pMg, I, T):
        v_dG0 = []
        v_nH  = []
        v_z   = []
        v_nMg  = []
        
        for key, dG0_list in self.dgs.iteritems():
            nH, z, nMg = key
            for dG0 in dG0_list:
                v_dG0.append(dG0)
                v_nH.append(nH)
                v_z.append(z)
                v_nMg.append(nMg)

        v_dG0 = pylab.array(v_dG0)
        v_nH  = pylab.array(v_nH)
        v_z   = pylab.array(v_z)
        v_nMg  = pylab.array(v_nMg)
        
        v_ddG0 = correction_function(v_nH, v_nMg, v_z, pH, pMg, I, T)
        v_dG0_tag = v_dG0 + v_ddG0

        return v_dG0, v_dG0_tag, v_nH, v_z, v_nMg 

    def Transform(self, pH, pMg, I, T):
        """
            Transform this set of pseudoisomers to a dG value at the 
            specified conditions.
        """
        if not self.dgs:
            raise ValueError("Cannot run Transform on an empty pseudoisomer map")

        _v_dG0, v_dG0_tag, _v_nH, _v_z, _v_nMg = self._Transform(pH, pMg, I, T)
        return -R * T * log_sum_exp(v_dG0_tag / (-R*T))
    
    def GetMostAbundantPseudoisomer(self, pH, I, pMg, T):
        v_dG0, v_dG0_tag, v_nH, v_z, v_nMg = self._Transform(pH, pMg, I, T)
        i = pylab.argmin(v_dG0_tag)
        return v_dG0[i], v_dG0_tag[i], v_nH[i], v_z[i], v_nMg[i]

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
        for (nH, z, nMg), dG0_list in self.dgs.iteritems():
            for dG0 in dG0_list:
                res.append((nH, z, nMg, dG0))
        return sorted(res)
    
    def __str__(self):
        s = ""
        for (nH, z, nMg), dG0_list in self.dgs.iteritems():
            s_dG0 = ", ".join(["%.1f" % dG0 for dG0 in dG0_list])
            s += "nH=%d z=%d nMg=%d dG0=%s kJ/mol\n" % (nH, z, nMg, s_dG0)
        return s

    def Display(self, cid):
        for nH, z, nMg, dG0 in self.ToMatrix():
            print "C%05d | %2d | %2d | %3d | %6.2f" % (cid, nH, z, nMg, dG0)
    
    def GetRef(self, nH, z, nMg):
        key = self._MakeKey(nH, z, nMg)
        return self.refs.get(key, None)

    def GetdG0(self, nH, z, nMg, T=default_T):
        key = self._MakeKey(nH, z, nMg)
        if key not in self.dgs:
            return None
        
        dG0_list = self.dgs[key]
        return -(R*T) * log_sum_exp([x/(-R*T) for x in dG0_list])
    
    def GetpKa(self, nH, z, nMg, T=default_T):
        dG0_f_deprotonated = self.GetdG0(nH-1, z-1, nMg, T)
        dG0_f_protonated = self.GetdG0(nH, z, nMg, T)
        if not dG0_f_deprotonated or not dG0_f_protonated:
            return None
        else:
            return (dG0_f_deprotonated - dG0_f_protonated) / (R*T*pylab.log(10))
        
    def GetAllpKas(self, nMg=0, T=default_T):
        """
            Returns:
                a list of tuples with (nH_below, nH_above, pKa) for each 
                protonation which could be calculated
        """
        res = []
        for curr_nH, curr_z, curr_nMg in self.dgs.keys():
            if nMg == curr_nMg:
                pKa = self.GetpKa(curr_nH, curr_z, curr_nMg, T)
                if pKa is not None:
                    res.append((curr_nH, curr_nH-1, pKa))
        return sorted(res, reverse=True)
        
    def GetpK_Mg(self, nH, z, nMg, T=default_T):
        dG0_f_without_Mg = self.GetdG0(nH, z-2, nMg-1, T)
        dG0_f_with_Mg = self.GetdG0(nH, z, nMg, T)
        if not dG0_f_without_Mg or not dG0_f_with_Mg:
            return None
        else:
            return (dG0_f_without_Mg + dG0_f_Mg - dG0_f_with_Mg) / (R*T*pylab.log(10))
