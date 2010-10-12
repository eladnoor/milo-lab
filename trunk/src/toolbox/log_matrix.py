from pylab import log, exp, zeros, ones, matrix, pi, nan
from toolbox.util import log_sum_exp, log_subt_exp

def log_regress(X, y):
    return log_dot(log_dot(log_inv(log_dot(X.T,X)), X.T), y)

def log_mat(X):
    res = zeros(X.shape) * complex(0,0)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            if (X[i,j] > 0):
                res[i,j] = complex(log(X[i,j]), 0)
            elif (X[i,j] < 0):
                res[i,j] = complex(log(-X[i,j]), pi)
            else:
                res[i,j] = nan
    return res

def log_dot(X,Y):
    if (X.shape[1] != Y.shape[0]):
        raise Exception("The dimensions of the matrices do not match: (%d,%d) and (%d,%d)" % (X.shape[0], X.shape[1], Y.shape[0], Y.shape[1]))
    R = zeros((X.shape[0], Y.shape[1]))
    for i in range(X.shape[0]):
        for j in range(Y.shape[1]):
            R[i,j] = log_sum_exp([(X[i,k] + Y[k,j]) for k in range(X.shape[1])])
    return R

def log_det(X):
    if (X.shape[0] != X.shape[1]):
        raise Exception("X is not a square matrix so its determinant cannot be calculated")
    
    if (X.shape[0] == 1):
        return X[0,0]
    
    if (X.shape[0] == 2): # X is a 2x2 matrix
        return(log_subt_exp(X[0,0]+X[1,1], X[1,0]+X[0,1]))
    
    if (X.shape[0] == 3):    
        s_plus = log_sum_exp([X[0,0]+X[1,1]+X[2,2], X[0,2]+X[1,0]+X[2,1], X[0,1]+X[1,2]+X[2,0]])
        s_minus = log_sum_exp([X[0,1]+X[1,0]+X[2,2], X[0,2]+X[1,1]+X[2,0], X[0,0]+X[1,2]+X[2,1]])
        return log_subt_exp(s_plus, s_minus)

    raise Exception("log_det is only implemented for matrices of size < 4")

def log_inv(X): # inverts a 3x3 matrix given by the logscale values
    if (X.shape[0] != X.shape[1]):
        raise Exception("X is not a square matrix and cannot be inverted")
    
    if (X.shape[0] == 1):
        return matrix((-X[0,0]))
    
    ldet = log_det(X)
    if (ldet == nan):
        raise Exception("The determinant of X is 0, cannot calculate the inverse")
     
    if (X.shape[0] == 2): # X is a 2x2 matrix
        I = (-log_det(X)) * ones((2,2))
        I[0,0] += X[1,1]
        I[0,1] += X[0,1] + complex(0, pi)
        I[1,0] += X[1,0] + complex(0, pi)
        I[1,1] += X[0,0]
        return I
    
    if (X.shape[0] == 3): # X is a 3x3 matrix
        I = (-log_det(X)) * ones((3,3))
        I[0,0] += log_subt_exp(X[1,1]+X[2,2], X[1,2]+X[2,1])
        I[0,1] += log_subt_exp(X[0,2]+X[2,1], X[0,1]+X[2,2])
        I[0,2] += log_subt_exp(X[0,1]+X[1,2], X[0,2]+X[1,1])
        I[1,0] += log_subt_exp(X[2,0]+X[1,2], X[1,0]+X[2,2])
        I[1,1] += log_subt_exp(X[0,0]+X[2,2], X[0,2]+X[2,0])
        I[1,2] += log_subt_exp(X[0,2]+X[1,0], X[0,0]+X[1,2])
        I[2,0] += log_subt_exp(X[1,0]+X[2,1], X[2,0]+X[1,1])
        I[2,1] += log_subt_exp(X[2,0]+X[0,1], X[0,0]+X[2,1])
        I[2,2] += log_subt_exp(X[0,0]+X[1,1], X[0,1]+X[1,0])
        return I
    
    raise Exception("log_inv is only implemented for matrices of size < 4")

if (__name__ == '__main__'):
    X = matrix([[2,1,-5],[3,1,2],[3,4,1]])*10
    M = log_dot(X.T,X)
    print M
    print log_inv(M)
    print exp(log_dot(M, log_inv(M)))
    
    
