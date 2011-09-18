import logging, os
import numpy as np
import matplotlib.pyplot as plt

from pygibbs.kegg import Kegg
from pygibbs.nist_regression import NistRegression
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg, default_T
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.sparse_kernel import SparseKernel
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.html_writer import HtmlWriter

def vector2string(v, cids, kegg):
    nonzero_columns = np.nonzero(abs(v) > 1e-10)[0]
    gv = " + ".join(["%g %s (C%05d)" % (v[j], kegg.cid2name(int(cids[j])), cids[j]) for j in nonzero_columns])
    return gv

def matrix2html(html_writer, A, cids, kegg):
    dict_list = []
    for i in xrange(A.shape[0]):
        dict_list.append({'row':'%d' % i, 'reaction':vector2string(A[i,:], cids, kegg)})
    html_writer.write_table(dict_list, headers=['row', 'reaction'])

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

    # Non-Alberty zeros:
    fixed_cids[101] = (21, 0.0) # THF - not from Alberty

    # Alberty zeros:
    fixed_cids[3]   = (26, 0.0) # NAD(ox)
    fixed_cids[10]  = (32, 0.0) # CoA
    fixed_cids[127] = (30, 0.0) # glutathione(ox)
    fixed_cids[376] = (28, 0.0) # retinal(ox)
    
    # Directly measured values
    fixed_cids[4]   = (27, 22.65) # NAD(red) -- relative to NAD(ox)
    fixed_cids[147] = (5, 313.40) # adenine
    #fixed_cids[121] = (10, -738.79) # D-ribose


    # Alberty zeros which are not in NIST:
    #fixed_cids[524] = ( 0, 0.0) # cytochrome c(ox)
    #fixed_cids[16]  = (31, 0.0) # FAD(ox)
    #fixed_cids[139] = ( 0, 0.0) # ferredoxin(ox)
    #fixed_cids[61]  = (19, 0.0) # FMN(ox)
    #fixed_cids[343] = ( 0, 0.0) # thioredoxin(ox)
    #fixed_cids[399] = (90, 0.0) # ubiquinone(ox)
    
    public_db = SqliteDatabase("../data/public_data.sqlite")
    alberty = PsuedoisomerTableThermodynamics.FromDatabase(public_db, 
        'alberty_pseudoisomers', label=None, name='Alberty')
    alberty_cid2dG0 = {}
    alberty_cid2nH = {}
    for cid in alberty.get_all_cids():
        pmap = alberty.cid2PseudoisomerMap(cid)
        dG0, _dG0_tag, nH, _z, _nMg = pmap.GetMostAbundantPseudoisomer(
            pH=default_pH,I=default_I, pMg=default_pMg, T=default_T)
        alberty_cid2nH[cid] = nH
        alberty_cid2dG0[cid] = dG0
    
    if not os.path.exists(prefix + 'S.txt'):
        db = SqliteDatabase("../res/gibbs.sqlite")
        nist_regression = NistRegression(db)
        
        cid2nH = {}
        for cid in nist_regression.nist.GetAllCids():
            if cid in fixed_cids:
                cid2nH[cid] = fixed_cids[cid][0]
            elif cid in alberty_cid2nH:
                cid2nH[cid] = alberty_cid2nH[cid]
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
    else:
        C = np.loadtxt(prefix + 'CID.txt', delimiter=',')
        cids = [int(cid) for cid in C[:,0]]
        cid2nH = {}
        for i, cid in enumerate(cids):
            cid2nH[cid] = int(C[i, 1])
        S = np.loadtxt(prefix + 'S.txt', delimiter=',')
        dG0 = np.loadtxt(prefix + 'dG0.txt', delimiter=',')

    html_writer = HtmlWriter('../res/regression_fast.html')
    html_writer.write("<h1>Pseudoisomeric Reactant Contributions</h1>\n")
    html_writer.write("<p>The stoichiometric matrix (S):")
    html_writer.insert_toggle(start_here=True)
    html_writer.write_ul(['%d rows' % S.shape[0], '%d columns' % S.shape[1],
                          '%d rank' % np.linalg.matrix_rank(S)])
    matrix2html(html_writer, S, cids, kegg)
    html_writer.div_end()
    html_writer.write('</p>')
    
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
    cid2dG0 = {}
    for i, cid in enumerate(cids):
        cid2dG0[cid] = x[i]
    
    # Calculate the Kernel of the reduced stoichiometric matrix (after removing 
    # the columns of the fixed compounds). 
    cids_red = [cid for cid in cids if cid not in fixed_cids]
    index_red = [i for i in xrange(len(cids)) if i not in index2value]
    S_red = S[:, index_red]
    K_red = LinearRegression.Kernel(S_red)
    
    #print "Reduced Stoichiometric Matrix:"
    #print matrix2string(S_red, cids_red, kegg)
    #print '-'*80
    
    # Find all CIDs that are completely determined and do not depend on any
    # free variable. In other words, all zeros columns in K2.
    dict_list = []
    
    determined_indices = np.where(np.sum(abs(K_red), 0) < 1e-10)[0] # all zero-columns in reducedK
    determined_cids = [cids_red[i] for i in determined_indices]
    plot_data = []
    for cid in cids:
        d = {'CID':'C%05d' % cid, 'Compound':kegg.cid2name(cid),
             'nH':'%d' % cid2nH[cid], 'dG0 (PRC)':'%.1f' % cid2dG0[cid]}
        if cid in alberty_cid2dG0:
            d['dG0 (Alberty)'] = '%.1f' % alberty_cid2dG0[cid]
            if cid not in fixed_cids:
                plot_data.append((alberty_cid2dG0[cid], cid2dG0[cid], kegg.cid2name(cid)))
        else:
            d['dG0 (Alberty)'] = ''
        
        if cid in fixed_cids:
            d['Depends on'] = 'anchored'
        elif cid in determined_cids:
            d['Depends on'] = 'fixed compounds'
        else:
            d['Depends on'] = 'kernel dimensions'
        
        dict_list.append(d)
        
    dict_list.sort(key=lambda(x):(x['Depends on'], x['CID']))
    html_writer.write("<p>Formation energies determined by the linear constraints:")
    html_writer.insert_toggle(start_here=True)
    html_writer.write('<font size="1">')
    html_writer.write_table(dict_list, headers=['Compound', 'CID', 'nH', 'dG0 (PRC)', 'dG0 (Alberty)', 'Depends on'])
    html_writer.write('</font>')
    html_writer.div_end()
    html_writer.write('</p>')

    # Plot a comparison between PRC and Alberty formation energies    
    fig = plt.figure(figsize=(8,8), dpi=80)
    plt.plot([x[0] for x in plot_data], [x[1] for x in plot_data],
             'b.', figure=fig)
    for x, y, name in plot_data:
        plt.text(x, y, name, fontsize=6)
    plt.xlabel('Alberty $\Delta_f G^\circ$')
    plt.ylabel('PRC $\Delta_f G^\circ$')
    html_writer.write("<p>Plot comparing PRC and Alberty results:")
    html_writer.insert_toggle(start_here=True)
    html_writer.embed_matplotlib_figure(fig)
    html_writer.div_end()
    html_writer.write("</p>")
    
    K_sparse = SparseKernel(S_red).Solve()
    html_writer.write("<p>The null-space of the reduced stoichiometric matrix:")
    html_writer.insert_toggle(start_here=True)
    matrix2html(html_writer, K_sparse, cids_red, kegg)
    html_writer.div_end()
    html_writer.write("</p>")

                
if __name__ == "__main__":
    main()
