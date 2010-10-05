from pylab import *

def solve_pH(pKa, Ca, pKb, Cb, pKw=14):
    Ka = 10**(-pKa)
    Kb = 10**(-pKb)
    Kw = 10**(-pKw) # water constant

    a_4 = Kb
    a_3 = Kb * Ka + Kw + Kb * Cb
    a_2 = Kw * Ka + Cb * Kb * Ka - Kb * Kw - Ca * Ka * Kb
    a_1 = -(Kw * Ka * Kb + Kw**2 + Ca * Kw * Ka)
    a_0 = -Kw**2 * Ka

    p = poly1d([a_4, a_3, a_2, a_1, a_0])
    pH_vec = []
    for r in p.r:
        if (r > 0):
            pH_vec.append(-log10(r))
    
    if (len(pH_vec) == 1):
        return pH_vec[0]
    elif (len(pH_vec) == 0):
        raise Exception("No solution was found for the pH")
    else:
        sys.stderr.write('WARNING: more than one solution, returning the one closest to pH = 7\n');
        i = argmin(abs([pH - 7 for pH in pH_vec]))
        return pH_vec[i]

if (__name__ == "__main__"):
    Ca = 10**(-14) # Acid buffer concentration in M
    pKa = 7 # pKa of Acid
    pKb = 10 # pKb of Base
    pKw = 14

    Cb_vec = [10**i for i in arange(-4, 4, 0.1)]
    H_vec = []
    OH_vec = []
    A_vec = []
    BH_vec = []
    for Cb in Cb_vec:
        pH = solve_pH(pKa, Ca, pKb, Cb)
        H_vec.append(10**(-pH))
        OH_vec.append(10**(pH - pKw))
        A_vec.append(Ca / (1 + 10**(pKa - pH)))
        BH_vec.append(Cb / (1 + 10**(pKb + pH - pKw)))

    hold(True)
    plot(Cb_vec, H_vec, 'b')
    plot(Cb_vec, OH_vec, 'b--')
    plot(Cb_vec, A_vec, 'r--')
    plot(Cb_vec, BH_vec, 'g')
    xscale('log')
    yscale('log')
    xlabel('Total base concentration [M]')
    xlabel('ion concentration [M]')
    legend(['[H+]', '[OH-]', '[A-]', '[BH+]'])
    show()
