import csv
from pylab import log, sqrt, array
from toolbox.util import log_sum_exp

R = 8.31e-3 # kJ/(K*mol)
J_per_cal = 4.184
default_T = 298.15 # K
default_I = 0.1 # mM
default_pH = 7.0
default_c0 = 1 # M

class MissingCompoundFormationEnergy(Exception):
    def __init__(self, value, cid=0):
        self.value = value
        self.cid = cid
    def __str__(self):
        return repr(self.value)    

class Thermodynamics:
    def __init__(self):
        self.pH = default_pH
        self.I = default_I
        self.T = default_T
        
        self.c_mid = 1e-4
        self.c_range = (1e-6, 1e-2)
        self.bounds = {}

    @staticmethod
    def debye_huckel(I):
        return (2.91482 * sqrt(I)) / (1 + 1.6 * sqrt(I))

    @staticmethod
    def correction_function(nH, z, pH, I, T):
        """
            nH and z - are the species parameters (can be vectors)
            pH and I - are the conditions, must be scalars
            returns the correction element used in the transform function
        """
        DH = Thermodynamics.debye_huckel(I) / (R*T)
        return -nH * (log(10)*pH + DH) + (z**2) * DH

    @staticmethod
    def transform(dG0, nH, z, pH, I, T):
        return dG0 - R*T*Thermodynamics.correction_function(nH, z, pH, I, T)

    @staticmethod
    def array_transform(dG0, nH, z, pH, I, T):
        """
            dG0, nH and z - are the species parameters (can be vectors)
            pH and I - are the conditions, must be scalars
            returns the transformed gibbs energy: dG0'
        """
        return -(R*T) * log_sum_exp(dG0 / (-R*T) + Thermodynamics.correction_function(nH, z, pH, I, T))
    
    @staticmethod
    def pmap_to_dG0(pmap, pH, I, T, most_abundant=False):
        if (len(pmap) == 0):
            raise Exception("Empty pmap")
        
        v_dG0 = []
        v_nH  = []
        v_z   = []
        for ((nH, z), dG0_list) in pmap.iteritems():
            for dG0 in dG0_list:
                v_dG0.append(dG0)
                v_nH.append(nH)
                v_z.append(z)

        v_dG0 = array(v_dG0)
        v_nH  = array(v_nH)
        v_z   = array(v_z)

        if (most_abundant):
            return min(v_dG0 / (-R*T) + Thermodynamics.correction_function(v_nH, v_z, pH, I, T))
        else:
            return Thermodynamics.array_transform(v_dG0, v_nH, v_z, pH, I, T)

    def cid2pmap(self, cid):
        raise Exception("method not implemented")

    def get_all_cids(self):
        raise Exception("method not implemented")
        
    def cid_to_dG0(self, cid, pH=None, I=None, T=None):
        pH = pH or self.pH
        I = I or self.I
        T = T or self.T
        return Thermodynamics.pmap_to_dG0(self.cid2pmap(cid), pH, I, T, most_abundant=False)
    
    def reaction_to_dG0(self, sparse_reaction, pH=None, I=None, T=None):
        """
            calculate the predicted dG0_r
        """
        return sum([coeff * self.cid_to_dG0(cid, pH, I, T) for (cid, coeff) in sparse_reaction.iteritems()])
    
    def cid_to_bounds(self, cid, use_default=True):
        (curr_c_min, curr_c_max) = self.bounds.get(cid, (None, None))
        if (curr_c_min == None and use_default):
            curr_c_min = self.c_range[0]
        if (curr_c_max == None and use_default):
            curr_c_max = self.c_range[1]
        return (curr_c_min, curr_c_max)
            
    @staticmethod
    def pmap_to_matrix(pmap):
        res = []
        for ((nH, z), dG0_list) in pmap.iteritems():
            for dG0 in dG0_list:
                res.append((nH, z, dG0))
        return res

    @staticmethod
    def pmap_to_table(pmap):
        s = ""
        s += "%2s | %2s | %s\n" % ("nH", "z", "dG0_f")
        s += "-" * 20 + "\n"
        for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(pmap):
            s += "%2d | %2d | %.1f\n" % (nH, z, dG0)
        return s     

    def display_pmap(self, cid):
        for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(self.cid2pmap(cid)):
            print "C%05d | %2d | %2d | %6.2f" % (cid, nH, z, dG0)
    
    def write_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'nH', 'charge', 'dG0'])
        for cid in self.get_all_cids():
            for ((nH, z), dG0_list) in self.cid2pmap(cid).iteritems():
                for dG0 in dG0_list:
                    writer.writerow([cid, nH, z, dG0])
                
    def write_transformed_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'pH', 'I', 'T', 'dG0_tag'])
        for cid in self.get_all_cids():
            dG0_tag = self.cid_to_dG0(cid)
            writer.writerow([cid, self.pH, self.I, self.T, dG0_tag])