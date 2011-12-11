from toolbox.util import log_sum_exp
import csv
import numpy as np

R = 8.31e-3 # kJ/(K*mol)
F = 96.48 # kC/mol
J_per_cal = 4.184
default_T = 298.15 # K
default_I = 0.25 # mM
default_pH = 7.0
default_c0 = 1 # M
default_pMg = 10
default_c_mid = 1e-3 # M
default_c_range = (1e-6, 1e-2) # M
dG0_f_Mg = -455.3 # kJ/mol, formation energy of Mg2+

def debye_huckel(I):
    return (2.91482 * np.sqrt(I)) / (1 + 1.6 * np.sqrt(I))

def correction_function(nH, z, nMg, pH, pMg, I, T):
    """
        nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the correction element used in the transform function
        
    Returns:
        The correction, in units of kJ/mol.
    """
    DH = debye_huckel(I)
    return nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + nH * (R*T*np.log(10)*pH + DH) - (z**2) * DH

def transform(dG0, nH, z, nMg, pH, pMg, I, T):
    return dG0 + correction_function(nH, z, nMg, pH, pMg, I, T)

def array_transform(dG0, nH, z, nMg, pH, pMg, I, T):
    """
        dG0, nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the transformed gibbs energy: dG0'
    """
    ddG0 = correction_function(nH, z, nMg, pH, pMg, I, T)
    dG0_tag = dG0 + ddG0
    return -R * T * log_sum_exp(dG0_tag / (-R*T))


class RedoxCarrier(object):

    # dG0 =  -E'*F * deltaE - R*T*ln(10)*pH * deltaH
    # Where: 
    #    F = 96.48 # kC/mol
    #    R*T*ln(10) = 5.7 kJ/mol
    #    deltaE - change in e-
    #    deltaH - change in H+
    #    pH - the conditions in which the E' was measured
    
    def __init__(self, cid_ox, cid_red, nH_ox, nH_red, z_ox, z_red, E_prime, pH, ref):
        self.cid_ox = cid_ox
        self.cid_red = cid_red
        self.nH_ox = nH_ox
        self.nH_red = nH_red
        self.z_ox = z_ox
        self.z_red = z_red
        self.E_prime = E_prime
        self.pH = pH
        self.ref = ref
        self.delta_H = nH_red - nH_ox
        self.delta_e = (nH_red - nH_ox) - (z_red - z_ox) # difference in no. of electrons
        self.ddG0 = -E_prime * F * self.delta_e - \
                     R * default_T * np.log(10) * pH * self.delta_H

class RedoxCarriers(dict):
    
    def __init__(self):
        for row in csv.DictReader(open('../data/thermodynamics/redox.csv', 'r')):
            name = row['name']
            cid_ox = int(row['CID_ox'])
            cid_red = int(row['CID_red'])
            nH_ox = int(row['nH_ox'])
            z_ox = int(row['charge_ox'])
            nH_red = int(row['nH_red'])
            z_red = int(row['charge_red'])
            E_prime = float(row['E_tag'])
            pH = float(row['pH'])
            ref = row['ref']
            self[name] = RedoxCarrier(cid_ox, cid_red, nH_ox, nH_red, 
                                      z_ox, z_red, E_prime, pH, ref)
