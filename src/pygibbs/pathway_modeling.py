#/usr/bin/python2.6

import cvxmod
from cvxmod import atoms
import math
import numpy

R = 8.31e-3 # kJ/(K*mol)
T = 298.15 # K
RT = R*T


class Pathway(object):
    
    def __init__(self, S, formation_energies):
        """Create a pathway object.
        
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            formation_energies: formation energies for the compounds
                in standard conditions.
        """
        self.S = S
        self.dG0_f = formation_energies
        self.Nr, self.Nc = S.shape
        
    def FindpCr(self, c_mid=1e-3, ratio=3.0,
                regularized=False):
        """Compute the pCr using hip-hoptimzation!"""
        # Make the formation energy variables and bounds.
        dgf_primes = cvxmod.optvar('dGf', self.Nc)
        formation_lb = cvxmod.matrix([-1e6]*self.Nc)
        formation_ub = cvxmod.matrix([1e6]*self.Nc)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr to flow forward.
        reaction_lb = cvxmod.matrix([-1e3]*self.Nr)
        reaction_ub = cvxmod.matrix([0]*self.Nr)
        
        # Compute formation energies at c_mid.
        to_mid = lambda x: x + R*T*math.log(c_mid)
        dg_mids = cvxmod.matrix(map(to_mid, self.dG0_f))
        r_frac_1 = ratio / (ratio + 1.0)
        r_frac_2 = 1.0 / (ratio + 1)
        ln10 = math.log(10)
        
        # Make the objective and problem.
        pc = cvxmod.optvar('pC', 1)
        obj = pc
        if regularized:
            obj += cvxmod.norm2(dgf_primes - dg_mids)
        problem = cvxmod.problem(objective=cvxmod.minimize(obj))
        
        # Set the constraints.
        S = cvxmod.matrix(self.S)
        problem.constr.append(dgf_primes >= formation_lb)
        problem.constr.append(dgf_primes <= formation_ub)
        problem.constr.append(S*dgf_primes <= reaction_ub)
        problem.constr.append(S*dgf_primes >= reaction_lb)
        problem.constr.append(dgf_primes + r_frac_1 * RT * ln10 * pc >= dg_mids)
        problem.constr.append(dgf_primes - r_frac_2 * RT * ln10 * pc <= dg_mids)
        #print problem

        problem.solve()
        
        print "Optimal problem value is %.4f." % problem.value
        cvxmod.printval(pc)
        opt_dgs = cvxmod.value(dgf_primes)
        print 'min dGf', min(opt_dgs)
        print 'max dGf', max(opt_dgs)
        print opt_dgs
        
if __name__ == '__main__':
    S = numpy.array([[-1,1,0],[0,-1,1]])
    cvxmod.randseed()
    dGs = cvxmod.randn(1,3,std=1000)
    print dGs
    path = Pathway(S, dGs)
    path.FindpCr()
    path.FindpCr(regularized=True)
