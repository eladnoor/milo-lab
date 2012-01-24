import numpy as np
import matplotlib.pyplot as plt

def simulate(S0, Sn, k_plus, k_minus):
    n = len(k_plus)
    assert len(k_minus) == n

    S = np.zeros((n+1))
    S[0] = S0
    S[n] = Sn
    
    # shift k's so the index will start at 1
    k_plus = np.array([1.0] + k_plus, dtype='float')
    k_minus = np.array([1.0] + k_minus, dtype='float')
    q = k_plus/k_minus
    Q = lambda i,j : np.prod(q[i:j+1])
    
    J = (S[0] * Q(1,n) - S[n]) / (np.sum([Q(i,n)/k_plus[i] for i in xrange(1, n+1)]))
    
    for i in xrange(1, n+1):
        S[i] = (S[i-1] - J/k_plus[i]) * q[i]

    return J, S

def calc(k_1, k_2, k_m2):
    S0, Sn = 1.0, 1.0
    Q = 100.0
    k_m1 = k_1*k_2/(k_m2 * Q)
    J, S = simulate(S0, Sn, [k_1, k_2], [k_m1, k_m2])
    return J, S[1]

if __name__ == "__main__":
    ks = [10**x for x in np.arange(-3, 3, 0.1)]
    
    Js = []
    S1s = []
    for k in ks:
        #J, S1 = calc(10, 1, k_m2)
        J, S1 = calc(k, 1, 0.1)
        Js.append(J)
        S1s.append(S1)

    fig = plt.figure()
    plt.plot(ks, Js, figure=fig, label='$J$')
    plt.plot(ks, S1s, figure=fig, label='$[S_1]$')
    plt.xscale('log', figure=fig)
    plt.yscale('log', figure=fig)
    #plt.xlabel('$k_{-2}$')
    plt.xlabel('$k_1$')
    plt.legend()
    plt.show()