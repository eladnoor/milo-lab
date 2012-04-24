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
default_RT = R * default_T
default_c_mid = 1e-3 # M
default_c_range = (1e-6, 1e-2) # M
dG0_f_Mg = -455.3 # kJ/mol, formation energy of Mg2+

symbol_d_G = "&Delta;G"
symbol_d_G0 = "&Delta;G&deg;"
symbol_d_G_prime = "&Delta;G'"
symbol_d_G0_prime = "&Delta;G'&deg;"

symbol_dr_G = "&Delta;<sub>r</sub>G"
symbol_dr_G0 = "&Delta;<sub>r</sub>G&deg;"
symbol_dr_G_prime = "&Delta;<sub>r</sub>G'"
symbol_dr_G0_prime = "&Delta;<sub>r</sub>G'&deg;"
symbol_dr_Gc_prime = "&Delta;<sub>r</sub>G'<sup>c</sup>"

symbol_df_G = "&Delta;<sub>f</sub>G"
symbol_df_G0 = "&Delta;<sub>f</sub>G&deg;"
symbol_df_G_prime = "&Delta;<sub>f</sub>G'"
symbol_df_G0_prime = "&Delta;<sub>f</sub>G'&deg;"

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
    from toolbox.util import log_sum_exp
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
    
    def __init__(self, cid_ox, cid_red, nH_ox, nH_red, z_ox, z_red, E0_prime, pH, ref):
        self.cid_ox = cid_ox
        self.cid_red = cid_red
        self.nH_ox = nH_ox
        self.nH_red = nH_red
        self.z_ox = z_ox
        self.z_red = z_red
        self.E0_prime = E0_prime
        self.pH = pH
        self.ref = ref
        self.delta_H = nH_red - nH_ox
        self.delta_e = (nH_red - nH_ox) - (z_red - z_ox) # difference in no. of electrons
        self.ddG0_prime = -E0_prime * F * self.delta_e
        
        # this calculation is not correct, one must use the reverse Lagendre transform
        # in order to convert G' to G.
        self.ddG0 = self.ddG0_prime - R * default_T * np.log(10) * pH * self.delta_H

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
            E0_prime = float(row["E'0"])
            pH = float(row['pH'])
            ref = row['ref']
            self[name] = RedoxCarrier(cid_ox, cid_red, nH_ox, nH_red, 
                                      z_ox, z_red, E0_prime, pH, ref)
