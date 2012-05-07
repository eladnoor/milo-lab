import numpy as np
import scipy


def ln_stirlings(n):
    return 0.5*np.log(2.0 * np.pi * n) + n * np.log(n) - n


def ln_stirling_binomial(n, k):
    return (ln_stirlings(n) - ln_stirlings(k)
            - ln_stirlings(n-k))


def CalcPValues(count_mat, prob_mat=None, exact=False):
    """Calculates p-values for the co-incidence matrix of two variables.
    
    Computes the probability we would observe a greater (i,j) value in mat
    by randomly sampling from the global distribution of row and column values.
    
    Suppose the rows are oxygen requirements (aerobe, anaerobe) and the
    columns are genotypes (ed, emp). Overall, there are 75 aerobes and 75 anaerobes,
    while there are 100 emp genotypes and 50 ed. If the genes and oxygen requirements
    were randomly distributed among organisms independently then the probability
    of randomly choosing an anaerobe with ed genes would be (50/150) * (75/150) =
    1/6 (meaning there should be 25 such organisms). However, suppose we observe
    that 1/4 of organisms (39) are anaerobes with ed genes. The probability of
    observing a more extreme value than 1/4 by drawing randomly is 
    
    sum_{i > 39} ((150 choose i) * 1/6^i * (5/6)^(150-i))
    
    which is the sum of the probability of randomly choosing any single value 
    greater than 39.
    
    Args:
        count_mat: matrix with independent variable values on rows, dependents on columns.
        prob_mat: the matrix with probabilities of each (i,j) point. If None, will
                  computed according to counts.
        exact: if exact binomial coefficients should be computed.
    
    Returns:
        A matrix with the same shape as mat with a p-value at each i,j position.
    """
    total = np.sum(count_mat)
    ceil_total = int(np.ceil(total))
    total2 = float(total**2)
    
    pval_mat = np.matrix(np.zeros(count_mat.shape))
    rows, cols = count_mat.shape
    for i in xrange(rows):
        for j in xrange(cols):
            # Get the total count
            observed_count = count_mat[i,j]
            
            # Compute or fetch the probability of this pair.
            if prob_mat is not None:
                prob = prob_mat[i,j]
            else:
                total_ind = np.sum(count_mat[i, :])
                total_dep = np.sum(count_mat[:, j])
                prob = float(total_ind * total_dep) / total2
            
            # Calculate a p-value for each more extreme case.
            pvals = []
            floor_count = int(np.floor(observed_count))
            for higher_count in xrange(floor_count+1, ceil_total+1):
                comb = scipy.comb(ceil_total, higher_count, exact=exact)
                # If we can't compute the actual value of the binomial coefficient
                # then approximate the log using stirlings approximation
                if not np.isfinite(comb):
                    comb = ln_stirling_binomial(ceil_total, higher_count)
                    
                pval = np.log(comb)
                pval += (higher_count*np.log(prob))
                pval += ((total-higher_count)*np.log(1-prob))
                pvals.append(np.exp(pval))
            # our p-value is the sum over all more extreme cases
            pval_mat[i,j] = np.sum(pvals)
    return pval_mat


def Test():
    # Double check stirlings approximation implementation
    test_vals = [(30, 5),
                 (18, 7),
                 (91, 32)]
    for n, k in test_vals:
        print np.log(scipy.comb(n, k))
        print ln_stirling_binomial(n, k)
    
    m = np.matrix([[50, 39],
                   [66, 5]])
    print CalcPValues(m)
    
    prob_m = np.matrix([[0.25, 0.25],
                        [0.25, 0.25]])
    print CalcPValues(m, prob_m)
    
    
if __name__ == '__main__':
    Test()