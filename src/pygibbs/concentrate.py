#!/usr/bin/python2.6
"""
    Calculate concentrations of compounds for given set of reactions
"""
import pylab
import getopt
import types
import glpk
import ctypes
import sys

def linprog(f, A, b, lb=None, ub=None, verbose=False):
    """
        The constraints are:
            A*x <= b
            Aeq*x = beq
            lb <= x <= ub
        
        The optimization rule is:
            minimize f*x
        
        In matrix A:
            rows are reactions (indexed by 'r') - total of Nr
            columns are compounds (indexed by 'c') - total of Nc
    """

    Nr = A.shape[0]
    Nc = A.shape[1]

    lp = glpk.glp_create_prob()
    glpk.glp_set_prob_name(lp, "bottlenecks")
    glpk.glp_set_obj_dir(lp, glpk.GLP_MIN)
    glpk.glp_add_rows(lp, Nr)
    
    for r in range(Nr):
        glpk.glp_set_row_name(lp, r+1, "Reaction %d" % r)
        glpk.glp_set_row_bnds(lp, r+1, glpk.GLP_UP, 0.0, b[r]) # 0.0 is ignored since the lower bound is -infty

    # Create the columns, name the reactions (RID) and add the values to a temporary sparse matrix
    glpk.glp_add_cols(lp, Nc)
    for c in range(Nc):
        glpk.glp_set_col_name(lp, c+1, "Compound %d" % c)
        if (lb != None and lb[c] != None and ub != None and ub[c] != None):
            glpk.glp_set_col_bnds(lp, c+1, glpk.GLP_DB, lb[c], ub[c])
        elif (lb != None and lb[c] != None):
            glpk.glp_set_col_bnds(lp, c+1, glpk.GLP_DN, lb[c], 0.0)
        elif (ub != None and ub[c] != None):
            glpk.glp_set_col_bnds(lp, c+1, glpk.GLP_UP, 0.0, ub[c])

    # Copy the sparse matrix to C-type structures that GLPK can later use
    size = Nr*Nc 
    ia = (ctypes.c_int    * (1+size))() # array of row indices in sparse matrix
    ja = (ctypes.c_int    * (1+size))() # array of col indices in sparse matrix
    ar = (ctypes.c_double * (1+size))() # array of values in sparse matrix
    for r in range(Nr):
        for c in range(Nc):
            ia[1+r*Nc+c] = ctypes.c_int(r+1)
            ja[1+r*Nc+c] = ctypes.c_int(c+1)
            ar[1+r*Nc+c] = ctypes.c_double(A[r,c])

    glpk.glp_load_matrix(lp, size, ia, ja, ar)
    
    for c in range(Nc):
        glpk.glp_set_obj_coef(lp, c+1, f[c])

    parm = glpk.glp_smcp()
    glpk.glp_init_smcp(parm)
    if (verbose):
        parm.msg_lev = glpk.GLP_MSG_ON
    else:
        parm.msg_lev = glpk.GLP_MSG_OFF
    retval = glpk.glp_simplex(lp, parm)
    if (retval != 0):
        return None
    
    f_min = glpk.glp_get_obj_val(lp)

    solution = pylab.zeros((Nc, 1))
    for c in range(Nc):
        solution[c,0] = glpk.glp_get_col_prim(lp, c+1)
    
    m = pylab.dot(A,solution) - b
    for i in range(Nr):
        if (m[i,0] > 1e-10):
            return None
    return solution

# stoichiometric matrix
S = pylab.array([\
[-1 ,-1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,0 ,-1  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0 ,-1  ,0  ,1 ,-1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0 ,-1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0  ,0 ,-1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0  ,0  ,0 ,-1 ,-1 ,-1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,1  ,0 ,-1  ,0  ,0  ,0  ,0  ,0  ,0 ,-1  ,0  ,1  ,0  ,0  ,0  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0 ,-1  ,1  ,0  ,0  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0 ,-1  ,1  ,0  ,0  ,1],\
[ 0  ,1  ,0 ,-1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0 ,-1  ,1  ,0  ,0],\
[ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,0  ,0 ,-1  ,0  ,0  ,0 ,-1  ,1  ,0]])

# deltaG0'-formation values for all compounds (in kJ/mol)
dG_f = pylab.matrix([-432.6,-2291.9,-1315.9,-1424.4,-1316.3,-2199.0,-1093.6,-1085.5,\
        1045.3,-1058.2,-2208.5,1108.2,-1352.9,-1352.9,-1185.1,-352.5,\
        -315.7,-156.8]).T

Nr = S.shape[0]
Nc = S.shape[1]
R = 8.31e-3 # gas constant (kJ/(K*mol))
T = 300 # temperature (K)
C_min = 1e-6 # lower bound on concentration
C_max = 1e-2 # upper bound on concentration
C_avg = 1e-3 # the average concentration

# compute right hand-side matrix,
# i.e. the deltaG0' of the reactions divided by -RT (unitless)
dG_r0 = pylab.dot(S, dG_f)

# compute the deltaG_m (i.e. taking into account the number of substrates and products)
x0 = pylab.ones((Nc, 1)) * pylab.log(C_avg) # the trivial assignment of concentrations
dG_rm = dG_r0 + R*T*pylab.dot(S, x0)

# function for minimization (x'*f)
f = pylab.ones((Nc, 1))
b = dG_r0/(-R*T)
initial_c = pylab.log(1e-4) # 0.1 mM

# now try using the new one-step system:
# first, add two new variables (columns) for the lower and upper margin
# and add 2*Nc rows for the concentration conditions (an upper and lower bound for each compound)
S1 = pylab.zeros((Nr+Nc*2, Nc+2))
b1 = pylab.zeros((Nr+Nc*2, 1))
S1[0:Nr,0:Nc] = S
b1[0:Nr,:] = b

# the goal vector will be to minimize the sum of the two margins (i.e. the allowed concentration gap)
f1 = pylab.zeros((Nc+2, 1))
f1[Nc,0] = 1
f1[Nc+1,0] = 1

ub1 = pylab.ones((Nc+2, 1)) * 0.0
ub1[Nc,0] = 100.0 # a very large number
ub1[Nc+1,0] = 100.0 # a very large number
lb1 = pylab.ones((Nc+2, 1)) * -100.0 # a very small number
lb1[Nc,0] = 0.0 # each margin is a positive number
lb1[Nc+1,0] = 0.0 # each margin is a positive number
for i in range(Nc):
    # add the constraint:    -log(c_i) - M_low < -log(initial_c)
    # i.e.   c_i > initial_c/M_low
    S1[Nr+2*i,   i]    = -1
    S1[Nr+2*i,   Nc]   = -1
    b1[Nr+2*i,   0]    = -initial_c
    # add the constraint:    log(c_i) - M_high < log(initial_c)
    # i.e.   c_i < initial_c*M_high
    S1[Nr+2*i+1, i]    =  1
    S1[Nr+2*i+1, Nc+1] = -1
    b1[Nr+2*i+1, 0]    =  initial_c

x1 = linprog(f1, S1, b1, lb1, ub1, verbose=False)
if (x1 == None):
    raise Exception("no solution was found")

dG_r = dG_r0 + R*T*pylab.dot(S, x1[0:Nc])
concentrations = pylab.exp(x1[0:Nc,0])

print "Concentrations (in mM): \n", ", ".join(["%.3f" % (c*1000) for c in concentrations])
print "dG_r[0] (in kJ/mol): \n", ", ".join(["%.2f" % (g/1000) for g in dG_r0])
print "dG_r[m] (in kJ/mol): \n", ", ".join(["%.2f" % (g/1000) for g in dG_rm])
print "dG_r[opt] (in kJ/mol): \n", ", ".join(["%.2f" % (g/1000) for g in dG_r])

pylab.figure()
pylab.ylabel('dG [kJ/mol]')
pylab.plot(pylab.cumsum(pylab.array(dG_r0) / 1000))
pylab.plot(pylab.cumsum(pylab.array(dG_rm) / 1000))
pylab.plot(pylab.cumsum(pylab.array(dG_r) / 1000))
pylab.legend(['dG\'_r (C=1M)', 'dG_r (C=1mM))', 'dG_r (C=opt)'])
pylab.show()
