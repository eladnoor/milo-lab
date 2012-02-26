#!/usr/bin/python

import pylab
import types
from thermodynamic_constants import R, default_T, correction_function
from toolbox.util import log_sum_exp
from pygibbs.thermodynamic_constants import dG0_f_Mg, default_I, default_pMg

class PseudoisomerMap(object):
    """A map from pseudoisomers to dG values."""
    
    def __init__(self, nH=None, z=None, nMg=None, dG0=None, ref=''):
        self.dgs = {}
        self.refs = {}
        
        if None not in [nH, z, nMg, dG0]:
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
        
    def FilterImprobablePseudoisomers(self, threshold=0.001, T=default_T):
        pH_range = pylab.arange(0, 14.01, 2)
        I_range = pylab.arange(0.0, 0.51, 0.1)
        pMg_range = pylab.arange(0, 14.01, 2)
        
        # iterate all physiological conditions and keep only psuedoisomers that
        # are relatively abundant in at least one condition
        species_to_keep = set()
        for pH in pH_range:
            for I in I_range:
                for pMg in pMg_range:
                    species2abundance = self._Bolzmann(pH, pMg, I, T)
                    species_to_keep.update([key for key, abundance in 
                                            species2abundance.iteritems() 
                                            if abundance > threshold])
        
        for key in set(self.dgs.keys()).difference(species_to_keep):
            del self.dgs[key]
            del self.refs[key]

    def TitrationPlot(self, I=default_I, pMg=default_pMg, T=default_T):
        pH_range = pylab.arange(0, 14.01, 0.2)
        
        # iterate all physiological conditions and keep only psuedoisomers that
        # are relatively abundant in at least one condition
        species_to_curve = {}
        for pH in pH_range:
            species2abundance = self._Bolzmann(pH, pMg, I, T)
            for key, abundance in species2abundance.iteritems():
                species_to_curve.setdefault(key, [])
                species_to_curve[key].append(abundance)

        fig = pylab.figure()
        for key, curve in species_to_curve.iteritems():
            pylab.plot(pH_range, curve, '-', label='nH=%d, z=%d, nMg=%d' % key,
                       figure=fig)
        pylab.legend()
        pylab.show()
        
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
        
        v_ddG0 = correction_function(v_nH, v_z, v_nMg, pH, pMg, I, T)
        v_dG0_tag = v_dG0 + v_ddG0

        return v_dG0, v_dG0_tag, v_nH, v_z, v_nMg
    
    def _Bolzmann(self, pH, pMg, I, T):
        """
            Returns a dictionary with the same keys as self.dgs, but the 
            values are the relative abundance of that species at the 
            specified conditions
        """
        _v_dG0, v_dG0_tag, v_nH, v_z, v_nMg = self._Transform(pH, pMg, I, T)
        minus_E_over_RT = [-E / (R*T) for E in v_dG0_tag]
        total = log_sum_exp(minus_E_over_RT)
        
        species2abundance = {}
        for i in xrange(len(v_dG0_tag)):
            key = PseudoisomerMap._MakeKey(v_nH[i], v_z[i], v_nMg[i])
            abundance = pylab.exp(minus_E_over_RT[i] - total)
            species2abundance[key] = species2abundance.get(key, 0.0) + abundance
        
        return species2abundance

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
        for key, dG0_list in self.dgs.iteritems():
            (nH, z, nMg) = key
            for dG0 in dG0_list:
                res.append((nH, z, nMg, dG0))
        return sorted(res)
    
    def WriteToHTML(self, html_writer):
        dict_list = []
        for key, dG0_list in self.dgs.iteritems():
            (nH, z, nMg) = key
            d = {'nH':nH, 'charge':z, 'nMg':nMg}
            d['dG0 [kJ/mol]'] = ', '.join(['%.1f' % dG0 for dG0 in dG0_list])
            d['reference'] = self.refs.get(key, '')
            dict_list.append(d)
        dict_list.sort(key=lambda(k):(k['nH'], k['nMg']))
        html_writer.write_table(dict_list, headers=['nH', 'charge', 'nMg',
                                                    'dG0 [kJ/mol]', 'reference'])        
        
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
    
    def SetRef(self, nH, z, nMg, ref):
        key = self._MakeKey(nH, z, nMg)
        self.refs[key] = ref

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


if __name__ == "__main__":
    pmap = PseudoisomerMap()
    pmap.Add(1, 0, 0, -150.0)
    pmap.Add(2, 1, 0, -200.0)
    pmap.Add(3, 2, 0, -210.0)
    
    print pmap
    print '*'*50
    print pmap.TitrationPlot()
    pmap.FilterImprobablePseudoisomers()
    print pmap
    print '*'*50
    