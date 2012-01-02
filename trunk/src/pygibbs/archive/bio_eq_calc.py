from pylab import matrix, array, dot, sqrt, zeros, ones, exp, log, inv
import scipy.linalg
import math, sys
from log_matrix import *

R = 8.31e-3 # kJ/mol
T = 300.0 # K

def debye_huckel(I):
    return (2.91482 * sqrt(I)) / (1 + 1.6 * sqrt(I))

def transform(species, pH, I):
    dG_tag = []
    for (nH, z, dG) in species:
        dG_tag.append(dG + R*T*nH*log(10)*pH - (z**2 - nH)*debye_huckel(I))
    dG_tag_total = -(R*T) * log_sum_exp(array(dG_tag) / (-R*T))
    return dG_tag_total

# no. of hydrogens (nH), charge (z), formation energy of species (in kJ/mol)
#species = [(12, -4, -2768.1), (13, -3, -2811.48), (14, -2, -2838.18)]
species = [(12, -4, -2768.1), (13, -3, -2811.48), (14, -2, -2838.18)]

conditions = [(4.0, 0.0), (5.0, 0.0), (7.0, 0.0), (8.0, 0.0), (6.0, 0.0)]
measurements = []
for (pH, I) in conditions:
    dG_tag = transform(species, pH, I)
    measurements.append((pH, I, dG_tag))

print "dG (original):\n", species
print "dG' (measurements):\n", measurements

Nm = len(measurements)
Ns = len(species)

# These are the definitions for 'y', 'X' and 'a':
#    y[i] = exp(-dG'[i] / RT)
#    X[i,j] = exp(-nH*pH*log(10)) + (z^2 - nH)*debye_huckel(I) / RT
#    a[j] = exp(-dG[i] / RT)
#    y = X*a

Ly = zeros((Nm,1))
LX = zeros((Nm, Ns))
for i in range(Nm):
    (pH, I, dG) = measurements[i]
    Ly[i] = -dG / (R*T)
    for j in range(Ns):
        (nH, z, dG) = species[j]
        LX[i,j] = -nH*pH*log(10) + (z**2 - nH)*debye_huckel(I) / (R*T)
    
    pivot = max(LX[i,:])
    LX[i,:] -= pivot
    Ly[i] -= pivot   


print "LX:\n", LX
print "Ly:\n", Ly
#LP_i = log_inv(LX)
#print "LP (inverse):\n", LP_i

Lcov = log_dot(LX.T,LX)
print "Lcov\n", Lcov
print "I\n", exp(log_dot(log_inv(Lcov),Lcov))

LP_r = log_dot(log_inv(Lcov), LX.T)
print "LP (regression):\n", LP_r

#La_i = log_dot(LP_i, Ly)
La_r = log_dot(LP_r, Ly)
#print    measurements.append((pH, I, dG_tag)) "La (inverse):\n", La_i
print "La (regression):\n", La_r

species_r = []
print "\n nH  |  z   | dG original | dG regress:"
for j in range(Ns):
    (nH, z, dG) = species[j]
    dG_r = -R*T*La_r[j,0]
    species_r.append((nH, z, dG_r)) 
    print "%4.1f | %4.1f | %11.1f | %11.1f" % (nH, z, dG, dG_r)

print "\n  pH |   I  | dG regress  | dG original:"
for i in range(Nm):
    (pH, I, dG_tag) = measurements[i]
    dG_tag_r = transform(species_r, pH, I)
    print "%4.1f | %4.1f | %11.1f | %11.1f" % (pH, I, dG_tag_r, dG_tag)
