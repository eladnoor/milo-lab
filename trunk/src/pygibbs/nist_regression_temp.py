import logging, os, sys
import numpy as np
from matplotlib import mlab

from pygibbs.kegg import Kegg
from pygibbs.nist_regression import NistRegression
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg, default_T
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.sparse_kernel import SparseKernel

def vector2string(v, cids, kegg):
    nonzero_columns = np.nonzero(abs(v) > 1e-10)[0]
    gv = " + ".join(["%g %s (C%05d)" % (v[j], kegg.cid2name(int(cids[j])), cids[j]) for j in nonzero_columns])
    return gv

def main():
    kegg = Kegg.getInstance()
    prefix = '../res/prc_'
    
    fixed_cids = {} # a dictionary from CID to pairs of (nH, dG0)
    
    # Alberty formation energies directly measured, linearly independent:
    fixed_cids[1]   = (2, -237.19) # H2O
    fixed_cids[9]   = (1, -1096.1) # HPO3(-2)
    fixed_cids[14]  = (4, -79.31) # NH4(+1)
    fixed_cids[59]  = (0, -744.53) # SO4(-2)
    fixed_cids[288] = (1, -586.77) # HCO3(-1)
    fixed_cids[147] = (5, 313.40) # adenine

    # Non-Alberty zeros:
    fixed_cids[101] = (21, 0.0) # THF - not from Alberty

    # Alberty zeros:
    fixed_cids[3]   = (26, 0.0) # NAD(ox)
    fixed_cids[10]  = (32, 0.0) # CoA
    fixed_cids[51]  = (16, 0.0) # glutathione(ox)
    fixed_cids[376] = (28, 0.0) # retinal(ox)

    # Alberty zeros which are not in NIST:
    #fixed_cids[524] = ( 0, 0.0) # cytochrome c(ox)
    #fixed_cids[16]  = (31, 0.0) # FAD(ox)
    #fixed_cids[139] = ( 0, 0.0) # ferredoxin(ox)
    #fixed_cids[61]  = (19, 0.0) # FMN(ox)
    #fixed_cids[343] = ( 0, 0.0) # thioredoxin(ox)
    #fixed_cids[399] = (90, 0.0) # ubiquinone(ox)
    
    if not os.path.exists(prefix + 'S.txt'):
        db = SqliteDatabase("../res/gibbs.sqlite")
        nist_regression = NistRegression(db)
        
        cid2nH = {}
        for cid in nist_regression.nist.GetAllCids():
            if cid in fixed_cids:
                cid2nH[cid] = fixed_cids[cid][0]
            else:
                tmp = nist_regression.dissociation.GetMostAbundantPseudoisomer(
                    cid, pH=default_pH, I=default_I, pMg=default_pMg, T=default_T)
                if tmp is not None:
                    cid2nH[cid] = tmp[0]
                else:
                    logging.warning('The most abundant pseudoisomer of %s (C%05d) '
                                    'cannot be resolved. Using nH = 0.' % (kegg.cid2name(cid), cid))
                    cid2nH[cid] = 0
        
        #nist_regression.std_diff_threshold = 2.0 # the threshold over which to print an analysis of a reaction
        #nist_regression.nist.T_range = None#(273.15 + 24, 273.15 + 40)
        S, dG0, cids = nist_regression.ReverseTransform(cid2nH=cid2nH)

        # export the raw data matrices to text files
        
        C = np.array([[cid, cid2nH.get(cid, 0)] for cid in cids])
        np.savetxt(prefix + 'CID.txt', C, fmt='%d', delimiter=',')
        np.savetxt(prefix + 'S.txt', S, fmt='%g', delimiter=',')
        np.savetxt(prefix + 'dG0.txt', dG0, fmt='%.2f', delimiter=',')
        
        for i in xrange(S.shape[0]):
            print i, NistRegression.ReactionVector2String(S[i, :], cids)

    else:
        C = np.loadtxt(prefix + 'CID.txt', delimiter=',')
        cids = [int(cid) for cid in C[:,0]]
        cid2nH = {}
        for i, cid in enumerate(cids):
            cid2nH[cid] = int(C[i, 1])
        S = np.loadtxt(prefix + 'S.txt', delimiter=',')
        dG0 = np.loadtxt(prefix + 'dG0.txt', delimiter=',')

    print "S has %d rows and %d columns and rank %d" % (S.shape[0], S.shape[1], np.linalg.matrix_rank(S))
    
    index2value = {}
    S_extended = S # the stoichiometric matrix, extended with elementary basis vector for the fixed compounds
    for cid in fixed_cids.keys():
        i = cids.index(cid)
        e_i = np.zeros((1, len(cids)))
        e_i[0, i] = 1.0
        S_extended = np.vstack([S_extended, e_i])
        nH, dG0_fixed = fixed_cids[cid]
        index2value[i] = dG0_fixed 
    x, K = LinearRegression.LeastSquaresWithFixedPoints(S, dG0, index2value)

    for cid in sorted(fixed_cids.keys()):
        i = cids.index(cid)
        nH, dG0_fixed = fixed_cids[cid]
        print "%40s (C%05d [nH=%2d]): dG0 = %.1f (should be %.1f)" % \
            (kegg.cid2name(cid), cid, cid2nH[cid], x[i], dG0_fixed)
    print '-'*80
    print "Null-space"

    # Find all CIDs that are completely determined and do not depend on any
    # free variable. In other words, all zeros columns in K2.

    K2 = SparseKernel(S_extended)
    underdetermined_indices = set()
    for i, v in enumerate(K2):
        nonzero_columns = np.where(abs(v) > 1e-10)[0]
        underdetermined_indices.update(nonzero_columns)
        print i, NistRegression.ReactionVector2String(v, cids)
    print '-'*80
    
    determined_cids = set(cids).difference([cids[i] for i in underdetermined_indices])
    for cid in sorted(determined_cids):
        i = cids.index(cid)
        print "%40s (C%05d [nH=%2d]): dG0 = %.1f" % \
            (kegg.cid2name(cid), cid, cid2nH[cid], x[i])
    
    sys.exit(0)
    
    # calculate the conservation matrix, i.e. a matrix with the known components
    # that are conserved throughout all the NIST reactions. The easiest are
    # the elements.
    elements = ['C', 'O', 'N', 'S', 'P']
    conv_mat = np.zeros((len(elements), S.shape[1]))
    for col, cid in enumerate(cids):
        atom_bag = kegg.cid2atom_bag(int(cid))
        if 'S' in atom_bag:
            print cid
        for row, atom in enumerate(elements):
            conv_mat[row, col] = atom_bag.get(atom, 0)

    conv_reactions = []            
    conv_reactions += [{3:1, 4:1, 5:1, 6:1}] # NAD(P)(H)
    conv_reactions += [{10:1, 24:1, 83:1, 91:1, 100:1, 154:1, 313:1, 332:1, 
                        683:1, 798:1, 920:1, 1144:1, 1213:1, 2232:1, 2557:1,
                        5268:1, 5269:1}] # CoA

    for conv_reaction in conv_reactions:
        conv_vector = np.zeros((1, S.shape[1]))
        for col, cid in enumerate(cids):
            conv_vector[0, col] = conv_reaction.get(cid, 0)
        conv_mat = np.vstack([conv_mat, conv_vector])
    
    for i in xrange(conv_mat.shape[0]):
        v = np.dot(S, conv_mat[i,:].T)
        if mlab.rms_flat(v) > 1e-6:
            print "ERROR: the %d conservation rule is not orthogonal to S" % i
            for j in xrange(S.shape[0]):
                if abs(v[j]) > 1e-6:
                    print j, vector2string(S[j, :], cids, kegg)
            sys.exit(-1)
        else:
            gv = vector2string(conv_mat[i,:], cids, kegg)
            print "conservation %d : %s" % (i, gv)
    
    # remove reactions that do not appear in any reaction
    nonzero_columns = np.sum(abs(S), 0).nonzero()[0]
    S = S[:, nonzero_columns]
    cids = [cids[i] for i in nonzero_columns]
        
    dG0 = dG0.reshape((dG0.shape[0], 1))
    
    K = SparseKernel(np.vstack([S, conv_mat]))
    #K = SparseKernel(S)
    
    for i, nullvector in enumerate(K):
        gv = vector2string(nullvector, cids, kegg)
        print "nullspace %d : %s" % (i, gv)
                
if __name__ == "__main__":
    main()
