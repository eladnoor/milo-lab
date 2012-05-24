
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
        self.inverse_kcat = 1 / self.kcat
        self.max_rate_factor = np.multiply(self.mass_per_mol_asite,
                                           self.inverse_kcat)
    
    def StoichFactor(self):
        return self.fluxes
    
    def MaximalRateFactor(self):
        return self.max_rate_factor

    def GetDGTag(self, full_x):
        return self.dG0 + RT * full_x * self.S
    
    def ThermoFactor(self, dgtag):
        return 1.0 / (1.0 - np.exp(dgtag/RT))
    
    def KineticFactor(self, full_x):
        exponent = np.log(self.km) - np.matrix(full_x).T
        exponent = np.multiply(self.m_plus, exponent)
        sum_exponents = np.sum(exponent, axis=0)
        exp_sum = np.exp(sum_exponents)
        return 1.0 + exp_sum

    def Derivatives(self, x):
        """Computes the derivative of our function.
        
        Our function f = sum(m_i u_i) = sum(C_u * T_u * K_u)
        where C is some constant, T is a thermodynamic factor 
        (a function of concentrations) and K is a kinetic factor
        (also function of concentrations).
        
        The derivative of u = (C * T * K) with respect to
        x_i = ln(c_i) is du/dx_i = C (T'K + K'T) where
           T' = n_i (T^2 - T)  and 
           K' = -m_i (K - 1)
        Note that n_i is the stoichiometric coefficient of metabolite
        i (negative for substrates, positive for products) and 
        m_i is the absolute stoichiometric coefficient if i is a 
        substrate and 0 otherwise.
        
        Args:
            x: the variable concentrations.
        
        Returns:
            A one-dimensional np.array of the Jacobian w.r.t. to x_i.
        """
        full_x = self.injector(x)
        
        dgtag = self.GetDGTag(full_x)
        T = self.ThermoFactor(dgtag)
        K = self.KineticFactor(full_x)

        # For numeric reasons, approximate where K 
        # is super large or when T is super small
        infinite_k = np.isinf(K)
        K[infinite_k] = 1e20
        small_T = np.where(np.abs(T) < 1e-20)
        T[small_T] = 0.0
        
        term_1 = np.multiply(-self.m_plus,
                             np.multiply((K-1.0), T))
        T_sq_min_T = np.square(T) - T
        term_2 = np.multiply(self.S,
                             np.multiply(T_sq_min_T, K))
        
        stoich_factor = self.StoichFactor()
        maximal_rate_factor = self.MaximalRateFactor()
        C = np.multiply(stoich_factor, maximal_rate_factor)
        
        du = np.multiply(C, (term_1 + term_2))
        return self.injector.Deject(np.sum(du, axis=1).T)
       
    def GetEnzymeLevels(self, x):
        """Get the modeled enzyme levels given concentrations."""
        my_x = self.injector(x)
        dgtag = self.GetDGTag(my_x)
        
        if (dgtag > 0).any():
            return np.matrix(np.ones(dgtag.shape) * 1e20)
        
        stoich_factor = self.StoichFactor()
        maximal_rate_factor = self.MaximalRateFactor()
        kinetic_factor = self.KineticFactor(my_x)
        thermo_factor = self.ThermoFactor(dgtag)
        
        try:
            const = np.multiply(stoich_factor, maximal_rate_factor)
            var = np.multiply(kinetic_factor, thermo_factor)
            levels = np.multiply(const, var)
        except Exception, e:
            print stoich_factor
            print maximal_rate_factor
            print kinetic_factor
            print thermo_factor
            raise e
        
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
    
    def Derivatives(self):
        raise NotImplementedError