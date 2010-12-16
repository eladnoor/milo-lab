from math import sqrt, log
from toolbox.util import log_sum_exp

R = 8.31e-3 # kJ/(K*mol)
J_per_cal = 4.184
default_T = 298.15 # K
default_I = 0.1 # mM
default_pH = 7.0
default_c0 = 1 # M
default_pMg = 3

dG0_f_Mg = -455.3 # kJ/mol, formation energy of Mg2+

def debye_huckel(I):
    return (2.91482 * sqrt(I)) / (1 + 1.6 * sqrt(I))

def correction_function(nH, nMg, z, pH, pMg, I, T):
    """
        nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the correction element used in the transform function
        
    Returns:
        The correction, in units of RT.
    """
    DH = debye_huckel(I) / (R*T)
    return nMg * (log(10)*pMg - dG0_f_Mg) + nH * (log(10)*pH + DH) - (z**2) * DH

def transform(dG0, nH, z, pH, I, T):
    return dG0 + R*T*correction_function(nH, z, pH, I, T)

def array_transform(dG0, nH, nMg, z, pH, pMg, I, T):
    """
        dG0, nH and z - are the species parameters (can be vectors)
        pH and I - are the conditions, must be scalars
        returns the transformed gibbs energy: dG0'
    """
    dG0_tag = dG0 + R*T*correction_function(nH, nMg, z, pH, pMg, I, T)
    return -R * T * log_sum_exp(dG0_tag / (-R*T))
