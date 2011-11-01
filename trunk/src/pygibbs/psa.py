import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from pygibbs.thermodynamic_constants import R, default_T

def phi(x, epsilon, kappa):
    if x <= 0:
        return np.inf
    return (np.sqrt(epsilon*(1+kappa)/x) + np.sqrt(epsilon*(1+kappa)/x + 1))**(-2)

def Phi(x, Epsilon, Kappa):
    """
        Input:
            x - float
            z - a list of floats (length = n)
        
        Calculate the function:
            p_i(x, z) = (sqrt(x/z[i]) + sqrt(x/z[i] + 1))^(-2)
            Phi(x, z) = prod(p_i)
            
    """
    return np.prod([phi(x, Epsilon[i], Kappa[i]) for i in xrange(len(Epsilon))])

def InvPhi(p, Epsilon, Kappa):
    """
        find x such that Phi(x, Z) = p
    """
    if p == 0:
        return 0.0
    if p < 0.01: # the values become very small when p -> 0, taking the reciprocal makes the search more stable
        roots = opt.fsolve(lambda(x):1.0/Phi(x, Epsilon, Kappa)-1.0/p, x0=p)
    else:
        roots = opt.fsolve(lambda(x):Phi(x, Epsilon, Kappa)-p, x0=p)
    return roots[0] # there should be only one, since Phi(x) is strongly monotonically decreasing

def EpsilonKappa(M, S, k_plus, k_minus):
    """
            M - list of enzyme Molecular Weights
            S - list of reaction stoichiometries (relative fluxes)
            k_plus - list of forward turnover numbers
            k_minus - list of backward turnover numbers
    """
    return ([M[i]*S[i] / k_plus[i]  for i in xrange(len(M))], 
            [k_plus[i] / k_minus[i] for i in xrange(len(M))])

def PSA1(Epsilon, verbose=False):
    """
        Use the method EpsilonKappa() to calculate Epsilon
    """
    if verbose:
        print "PSA1 costs:", ", ".join(["%.2g" % x for x in Epsilon])

    return 1.0 / np.sum(Epsilon)

def PSA2(Epsilon, Kappa, dG, T=default_T, verbose=False):
    """
        Use the method EpsilonKappa() to calculate Epsilon and Kappa
        
        dG - the overall change in Gibbs energy of the pathway
    """
    if dG == 0:
        return 0.0
    
    l = InvPhi(np.exp(dG/(R*T)), Epsilon, Kappa) # solve for the Lagrange multiplier
    
    costs = []
    for i in xrange(len(M)):
        costs.append(Epsilon[i] + l/(2 + 2*np.sqrt(1+l/(Epsilon[i]*(1+Kappa[i])))))
    
    if verbose:
        print "Epsilon: ", ", ".join(["%.2g" % x for x in Epsilon])
        print "Kappa:", ", ".join(["%.2g" % x for x in Kappa])
        print "PSA2 costs:", ", ".join(["%.2g" % x for x in costs])
        dG_steps = [R*T*np.log(phi(l, Epsilon[i], Kappa[i])) for i in xrange(len(Epsilon))]
        print "dGs:", ", ".join(["%.1f" % x for x in dG_steps])
        
    return 1.0 / np.sum(costs)

def PSA2_homogenous(N, epsilon, kappa, dG, T=default_T, verbose=False):
    phi = np.exp(dG/(N*R*T))
    
    costs = [epsilon * (1 + kappa*phi) / (1 - phi)] * N
    if verbose:
        print "Epsilon: ", ", ".join(["%.2g" % epsilon] * N)
        print "Kappa:", ", ".join(["%.2g" % kappa] * N)
        print "PSA2 costs:", ", ".join(["%.2g" % x for x in costs])
        print "dGs:", ", ".join(["%.1f" % phi] * N)

    return 1 / np.sum(costs)

if __name__ == "__main__":
    dG_vec = np.arange(-50, 1e-10, 2)

    fig = plt.figure(figsize=(8,8), dpi=90)
    plt.plot([min(dG_vec), max(dG_vec)], [1.0, 1.0], '--', linewidth=2, label='PSA$_1$')
    for m in [1e0, 1e1, 1e2, 1e3, 1e4]:    
        M = [0.5, 0.5]
        S = [1, 1]
        k_plus = [1.0, 1.0] # umol/min
        k_minus = [1.0/m, m] # umol/min
        
        Epsilon, Kappa = EpsilonKappa(M, S, k_plus, k_minus)
        
        psa = [PSA2(Epsilon, Kappa, dG) for dG in dG_vec] # in (umol/min)/mg
        plt.plot(dG_vec, psa, '-', linewidth=2, label='PSA$_2$ ($\kappa$ = %s)' % str(Kappa))
        
    plt.xlabel('pathway $\Delta G$')
    plt.ylabel('PSA')
    plt.legend(loc='lower left')
    plt.title('Pathway with 2 enzymes, both with a cost of $\epsilon = 0.5$')
    plt.ylim(0, 1.1)

    fig = plt.figure(figsize=(8,8), dpi=90)
    plt.plot([min(dG_vec), max(dG_vec)], [1.0, 1.0], '--', linewidth=2, label='PSA$_1$')
    for m in [1e2, 1e1, 1e0, 1e-1, 1e-2]:    
        M = [0.5, 0.5]
        S = [1, 1]
        k_plus = [1.0, 1.0] # umol/min
        k_minus = [m, m] # umol/min
        
        Epsilon, Kappa = EpsilonKappa(M, S, k_plus, k_minus)
        
        psa = [PSA2(Epsilon, Kappa, dG) for dG in dG_vec] # in (umol/min)/mg
        plt.plot(dG_vec, psa, '-', linewidth=2, label='PSA$_2$ ($\kappa$ = %s)' % str(Kappa))
        
    plt.xlabel('pathway $\Delta G$')
    plt.ylabel('PSA')
    plt.legend(loc='lower left')
    plt.title('Pathway with 2 enzymes, both with a cost of $\epsilon = 0.5$')
    plt.ylim(0, 1.1)
    plt.show()
