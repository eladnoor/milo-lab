import csv, re, random
from pygibbs.thermodynamics import Thermodynamics
from pygibbs.thermodynamic_constants import R, correction_function,\
    array_transform
from scipy.stats.stats import mean
from numpy.linalg.linalg import inv, LinAlgError, array
from matplotlib.pyplot import hold, figure, subplot, plot, xlabel, ylabel, title,\
    legend, show
from numpy.core.numeric import arange, zeros, exp
from toolbox.log_matrix import log_dot, log_mat
from toolbox.util import log_subt_exp

class ReverseTransformError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def solve(measurements, nH, z):
    """
        species - a list of pairs, containing the (nH, z) of the known species
        measurements - a list of 3-tuples, containing the (pH, I, delta-G) of the measurements
    """
    # Let 'X' be the auxiliary matrix:
    #    log(X[i,j]) = dG'[i] / RT - nH*pH*log(10)) + (z^2 - nH)*debye_huckel(I) / RT
    #
    # Let 'a' be the unknown column vector of the untransformed formation energies of the species:
    #    a[j] = exp(-dG[j] / RT)
    #
    # It is thus a fact that:
    #    1 = X*a
    #
    # So we shall use linear regression to find the values of a.

    Nm = len(measurements)
    Ns = len(nH)
    
    if (Nm == 0):
        raise Exception("There are no measurements given")
    if (Ns == 0):
        raise Exception("There are no species given")
  
    log_X = zeros((Nm, Ns))
    for i in range(Nm):
        (dG0_tag, pH, I, T) = measurements[i]
        log_X[i,:] = correction_function(nH, z, pH, I, T) + dG0_tag / (R*T)

    log_C = log_dot(log_X.T, log_X) # C = (X'*X)
    
    # To invert the matrix C, it is needed to normalize it first to be in a scale
    # where the inverse is numerically possible to calculate (dividing by the harmonic mean - M).
    # After inverting, we need to DIVIDE it again by the same factor, to return to the original scale.
    M = mean(log_C)
    try:
        log_inv_C = log_mat(inv(exp(log_C - M))) - M
    except LinAlgError:
        raise ReverseTransformError("Transform matrix is underdetermined (%d x %d), cannot solve pseudoisomer dG0" % (Nm, Ns))
    
    log_a = log_dot(log_inv_C, log_dot(log_X.T, zeros((Nm,1)))) # a = C_inv * X' * ones(Nm,1)
    
    dG0_solution = -R*T*array(log_a.flat)
    
    # This solution is usually biased (in linear-scale)
    # Therefore, we need to move its average in log-scale to fit the original data
    diff = zeros((Nm, 1))
    for i in range(Nm):
        (dG0_tag, pH, I, T) = measurements[i]
        dG_tag_pred = array_transform(dG0_solution, nH, z, pH, I, T)
        diff[i,0] = dG0_tag - dG_tag_pred
    
    return dG0_solution + mean(diff)

def test1():
    """
        Start from a set of known chemical formation energies (e.g. for ATP), calculate a set of random
        transformed dG0' values (add several levels of noise) and use the solver to recreate the chemical energies.
    """
    
    # no. of hydrogens (nH), charge (z), formation energy of species (in kJ/mol)
    nH = array([12, 13, 14])
    z = array([-4, -3, -2])
    dG0 = array([-2768.11, -2811.48, -2838.18])
    Ns = len(dG0)

    # conditions: pH, I, T
    #conditions = [(pH, 0.0, 300.0) for pH in [4.1, 4.5, 4.2, 4.3, 4.4, 4.5, 4.5, 6.5, 7.2, 8.9]]
    conditions = []
    for unused_i in range(1000):
        conditions.append((random.randrange(0.0, 18.0) + random.random(), 0.0, 300.0))
    
    figure()
    hold(True)
    plot_counter = 1
    for noise_level in [0.0, 0.5, 5.0, 20.0]: # in kJ/mol
        measurements = []
        noisy_measurements = []
        for (pH, I, T) in conditions:
            noise = random.normalvariate(0, noise_level)
            dG0_tag = array_transform(dG0, nH, z, pH, I, T)
            dG0_tag_noisy = dG0_tag+noise
            measurements.append((dG0_tag, pH, I, T))
            noisy_measurements.append((dG0_tag_noisy, pH, I, T))

            log_y = -dG0_tag / (R*T)
            log_y_noisy = -dG0_tag_noisy / (R*T)
            log_diff = log_subt_exp(log_y_noisy, log_y)
            print log_y, noise, log_diff
    
        dG0_r = solve(noisy_measurements, nH, z)
    
        print " nH  |  z   | dG regress  | dG error"
        print "--------------------------------------"
        for j in range(Ns):
            print "%4.1f | %4.1f | %11.2f | %11.2f" % (nH[j], z[j], dG0_r[j], dG0_r[j]-dG0[j])
        print "--------------------------------------"
        print ""
    
        #print ""
        #print "  pH |   I  | dG' original| dG' regress"
        #print "---------------------------------------"
        pH_list = arange(0.0, 18.001, 0.25)
        dG_tag_orig_list = []
        dG_tag_pred_list = []
        T = 300.0
        I = 0.0
        for pH in pH_list:
            dG_tag_orig = array_transform(dG0, nH, z, pH, I, T)
            dG_tag_pred = array_transform(dG0_r, nH, z, pH, I, T)
            #print "%4.1f | %4.1f | %11.2f | %11.2f" % (pH, I, dG_tag_orig, dG_tag_pred)
            dG_tag_orig_list.append(dG_tag_orig)
            dG_tag_pred_list.append(dG_tag_pred)
    
        subplot(2,2,plot_counter)
        plot([pH for (dG_noisy, pH, I, T) in noisy_measurements], [dG_noisy for (dG_noisy, pH, I, T) in noisy_measurements], 'r.')
        plot(pH_list, dG_tag_orig_list, 'g.')
        plot(pH_list, dG_tag_pred_list, 'b-')
        xlabel('pH')
        ylabel('dG')
        title('Extrapulated vs. Original (noise = %.1f)' % noise_level)
        legend(['Original', 'Measured (noisy)', 'Extrapulated'], loc=0)
        plot_counter += 1
        #plot(dG_tag_orig_list, dG_tag_orig_list, '-')
    show()

def test2():
    """
        Use Alberty's data (both chemical and biochemical tables) and do the reverse transform to try
        and recreate the chemical species formation energies"
    """
    chemical_file = open('../data/dG0_chemical.txt', 'r')
    compound2species = {}
    for line in chemical_file.readlines():
        try:
            [(compound_name, species_list)] = re.findall("([a-zA-Z0-9]+)sp=(.*)\n*", line)
        except ValueError as e:
            print str(e)
            raise e
        compound2species[compound_name] = eval(species_list)
    
    print "Reverse transform, according to changing pH\n" + "="*100

    csv_input = csv.reader(open("../data/dG0_biochemical_pH.csv", "r"))
    row = csv_input.next()
    pH = [5.0, 6.0, 7.0, 8.0, 9.0]
    I = 0.25
    T = 298.15
    for row in csv_input:
        compound_name = row[0]
        dG0 = [float(x) for x in row[1:]]
        measurements = []
        for i in range(5):
            measurements.append((dG0[i], pH[i], I, T))

        dG_array = array([dG for (dG, _, z, nH) in compound2species[compound_name]])
        z_array = array([z for (dG, _, z, nH) in compound2species[compound_name]])
        nH_array = array([nH for (dG, _, z, nH) in compound2species[compound_name]])
        
        dG0_r = solve(measurements, nH_array, z_array)
        print compound_name, ";".join(["(%.1f,%.1f)" % (dG_array[i], dG0_r[i]) for i in range(len(dG0_r))])

    print "\n\nReverse transform, according to changing Ionic strength\n" + "="*100

    csv_input = csv.reader(open("../data/dG0_biochemical_I.csv", "r"))
    row = csv_input.next()
    pH = 7.0
    I = [0.00, 0.10, 0.25]
    T = 298.15
    for row in csv_input:
        compound_name = row[0]
        dG0 = [float(x) for x in row[1:]]
        measurements = []
        for i in range(3):
            measurements.append((dG0[i], pH, I[i], T))

        dG_array = array([dG for (dG, _, z, nH) in compound2species[compound_name]])
        z_array = array([z for (dG, _, z, nH) in compound2species[compound_name]])
        nH_array = array([nH for (dG, _, z, nH) in compound2species[compound_name]])
        
        dG0_r = solve(measurements, nH_array, z_array)
        print compound_name, ";".join(["(%.1f,%.1f)" % (dG_array[i], dG0_r[i]) for i in range(len(dG0_r))])
    
def test3():
    T = 300
    dG0 = array([-2768.11, -2811.48, -2838.18])
    M = max(dG0)
    for unused_i in range(100):
        for j in range(dG0.shape[0]):
            e = random.normalvariate(0,1)
            log_y = -(dG0[j]-M)/(R*T)
            log_y_n = -(dG0[j]-M+e)/(R*T)
            diff = log_subt_exp(log_y_n, log_y)
            print diff
    
if (__name__ == '__main__'):
    test1()

