
import numpy as np

from pygibbs.thermodynamic_constants import default_T, R

RT = R * default_T

"""OLD VERSION
class EnzymeLevelFunc(object):
    
    def __init__(self, S, dG0, fluxes, kcat, km, masses,
                 injector):
        self.S = S
        self.m_plus = np.abs(np.clip(S, -1000, 0))
        self.dG0 = dG0
        self.fluxes = fluxes
        self.kcat = kcat
        self.km = km
        self.injector = injector
        self.scaled_fluxes = np.matrix(fluxes / kcat)
        self.Nr, self.Nc = self.S.shape
        self.mass_per_mol_asite = masses * 1000 # in g protein / mol
        #self.vtotal = 2 * 3.1e-6 # in mol s^-1 per GDW
        #self.vtotal = 8.2 / (60**2 * 1000) # in mol s^-1 per GDW
        self.vtotal = 5.0 / (60**2 * 1000) # in mol s^-1 per GDW
    
    def StoichFactor(self):
        return self.fluxes
    
    def MaximalRateFactor(self):
        inverse_kcat = 1 / self.kcat
        return np.multiply(self.mass_per_mol_asite, inverse_kcat)
    
    def ThermoFactor(self, full_x):
        dgtag = self.dG0 + RT * full_x * self.S
        if (dgtag >= 0).any():
            return np.matrix(np.ones(dgtag.shape) * 1e6) 
        return 1.0 / (1 - np.exp(dgtag/RT))
    
    def KineticFactor(self, full_x):
        x_exp = np.matrix(np.exp(full_x))
        scaled_kms = self.km / x_exp.T
        exponentiated = np.power(scaled_kms, self.m_plus)
        prods = np.prod(exponentiated, axis=0)
        return 1 + prods
    
    def GetEnzymeLevels(self, x):
        my_x = self.injector(x)
        
        stoich_factor = self.StoichFactor()
        maximal_rate_factor = self.MaximalRateFactor()
        kinetic_factor = self.KineticFactor(my_x)
        thermo_factor = self.ThermoFactor(my_x)
        
        levels = np.multiply(np.multiply(stoich_factor, maximal_rate_factor),
                             np.multiply(kinetic_factor, thermo_factor))
        #return levels
        return self.vtotal * levels * 100.0 * 100.0 / 55.0
        
    def __call__(self, x):
        levels = self.GetEnzymeLevels(x)
        return np.sum(levels)
"""

class ProteinCostFunc(object):
    """A functor computing optimal enzyme levels."""
    
    def __init__(self, S, dG0, fluxes, kcat, km, masses,
                 injector):
        self.S = S
        self.m_plus = np.abs(np.clip(S, -1000, 0))
        self.m_minus = np.abs(np.clip(S, 0, 1000))
        self.dG0 = dG0
        self.fluxes = fluxes
        self.kcat = kcat
        self.km = km
        self.injector = injector
        self.scaled_fluxes = np.matrix(fluxes / kcat)
        self.Nr, self.Nc = self.S.shape
        self.mass_per_mol_asite = masses * 1000 # in g protein / mol
    
    def StoichFactor(self):
        return self.fluxes
    
    def MaximalRateFactor(self):
        inverse_kcat = 1 / self.kcat
        return np.multiply(self.mass_per_mol_asite, inverse_kcat)
    
    def ThermoFactor(self, full_x):
        dgtag = self.dG0 + RT * full_x * self.S
        if (dgtag >= 0).any():
            return np.matrix(np.ones(dgtag.shape) * 1e6) 
        return 1.0 / (1 - np.exp(dgtag/RT))
    
    def KineticFactor(self, full_x):
        x_exp = np.matrix(np.exp(full_x))
        scaled_kms = self.km / x_exp.T
        exponentiated = np.power(scaled_kms, self.m_plus)
        prods = np.prod(exponentiated, axis=0)
        return 1 + prods
    
    def GetEnzymeLevels(self, x):
        """Get the modeled enzyme levels given concentrations."""
        my_x = self.injector(x)
        
        stoich_factor = self.StoichFactor()
        maximal_rate_factor = self.MaximalRateFactor()
        kinetic_factor = self.KineticFactor(my_x)
        thermo_factor = self.ThermoFactor(my_x)
        
        levels = np.multiply(np.multiply(stoich_factor, maximal_rate_factor),
                             np.multiply(kinetic_factor, thermo_factor))
        return levels
        
    def __call__(self, x):
        levels = self.GetEnzymeLevels(x)
        return np.sum(levels)
    
    

class ProteinCostFuncV2(ProteinCostFunc):
    """V2 functor models backwards reaction saturation."""
    
    def KineticFactor(self, full_x):
        x_exp = np.matrix(np.exp(full_x))
        scaled_concs = np.divide(x_exp.T, self.km)
        subs_exp = np.power(scaled_concs, self.m_plus)
        prods_exp = np.power(scaled_concs, self.m_minus)
        subs_prod = np.prod(subs_exp, axis=0)
        prods_prod = np.prod(prods_exp, axis=0)
        return np.divide((1.0 + subs_prod + prods_prod), subs_prod)