import numpy as np
from scipy.stats import pearsonr, spearmanr
import csv
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.kegg import Kegg 
from pygibbs.kegg_reaction import Reaction
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def LoadGrowthData():
    kegg = Kegg.getInstance()
    path = '../data/growth/growth_rates_adadi_2012.csv'
    data = []
    for row in csv.DictReader(open(path, 'r')):
        carbon_source = row['carbon source']
        cid, _, _ = kegg.name2cid(carbon_source)
        if cid is None:
            raise Exception("Cannot map compound name to KEGG ID: " + carbon_source)

        data.append({'carbon_source': carbon_source,
                     'cid': cid,
                     'growth_rate': float(row['maximum growth rate measured']),
                     'sumex': float(row['SUMEX'])})
    return data

def GetFullOxidationReaction(cid):
    kegg = Kegg.getInstance()

    basic_cids = [1, 7, 9, 11, 14] # H2O, O2, Pi, CO2, NH3
    basic_elements = ['C', 'O', 'P', 'N', 'e-']
    element_mat = np.matrix(np.zeros((len(basic_elements), len(basic_cids))))
    for j in xrange(len(basic_cids)):
        atom_bag = kegg.cid2atom_bag(basic_cids[j])
        atom_bag['e-'] = kegg.cid2num_electrons(basic_cids[j])
        for i in xrange(len(basic_elements)):
            element_mat[i, j] = atom_bag.get(basic_elements[i], 0)
            
    cs_element_vec = np.zeros((len(basic_elements), 1))
    atom_bag = kegg.cid2atom_bag(cid)
    atom_bag['e-'] = kegg.cid2num_electrons(cid)
    for i in xrange(len(basic_elements)):
        cs_element_vec[i, 0] = atom_bag.get(basic_elements[i], 0)
    
    x = np.linalg.inv(element_mat) * cs_element_vec
    
    sparse = dict([(basic_cids[i], np.round(x[i, 0], 3)) for i in xrange(len(basic_cids))])
    sparse[cid] = -1

    r = Reaction("complete oxidation of %s" % kegg.cid2name(cid), sparse)
    
    return r

def CorrPlot(x, y, labels, xlabel, ylabel, figure):
    plt.plot(x, y, '.', figure=figure)
    for j in xrange(len(x)):
        plt.text(x[j], y[j], labels[j], fontsize=8, figure=figure)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    (r_pearson, p_pearson) = pearsonr(x, y)
    (r_spearman, p_spearman) = spearmanr(x, y)
    plt.title(r'$r_{Pearson}$ = %.2f (p = %.3f), $r_{Spearman}$ = %.2f (p = %.3f)' %
              (r_pearson, p_pearson, r_spearman, p_spearman),
              fontsize=10)

def PlotEnergy(data, y_title, y_label, pdf):
    data_with_dG0 = [d for d in data if np.isfinite(d['dG0'])]

    plots = []
    plots.append({'x':[-d['dG0'] for d in data_with_dG0],
                  'y':[d[y_title] for d in data_with_dG0],
                  'labels':[d['carbon_source'] for d in data_with_dG0],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/mol]',
                  'ylabel':y_label,
                  })

    plots.append({'x':[-d['dG0']/d['numC'] for d in data_with_dG0],
                  'y':[d[y_title] for d in data_with_dG0],
                  'labels':[d['carbon_source'] for d in data_with_dG0],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/mol (C.S.)]',
                  'ylabel':y_label,
                  })

    plots.append({'x':[-d['dG0']/d['mw'] for d in data_with_dG0],
                  'y':[d[y_title] for d in data_with_dG0],
                  'labels':[d['carbon_source'] for d in data_with_dG0],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/gr]',
                  'ylabel':y_label,
                  })

    plots.append({'x':[d['total S'] for d in data],
                  'y':[d[y_title] for d in data],
                  'labels':[d['carbon_source'] for d in data],
                  'xlabel':r'$\Sigma_i |s_i|$',
                  'ylabel':y_label,
                  })
    
    fig = plt.figure(figsize=(12, 12), dpi=50)

    for i, d in enumerate(plots):
        plt.subplot(2,2,i+1)
        CorrPlot(d['x'], d['y'], d['labels'],
                 d['xlabel'], d['ylabel'], figure=fig)
    
    fig.tight_layout()
    
    pdf.savefig(fig)
    
def main():
    kegg = Kegg.getInstance()
    estimators = LoadAllEstimators()
    thermo = estimators['UGC']
    data = LoadGrowthData()
    
    pdf = PdfPages('../res/growth_rates.pdf')

    for d in data:
        mol = kegg.cid2mol(d['cid'])
        d['mw'] = mol.GetExactMass()
        atom_bag, _ = mol.GetAtomBagAndCharge()
        d['numC'] = atom_bag['C']
         
        r = GetFullOxidationReaction(d['cid'])
        r.Balance(balance_water=False, balance_hydrogens=True, exception_if_unknown=True)
        print r.FullReactionString(show_cids=False)
        d['dG0'] = r.PredictReactionEnergy(thermo)
        d['total S'] = sum([abs(x) for x in r.sparse.values()])
    
    PlotEnergy(data, 'growth_rate', 'Specific Growth Rate [1/hr]', pdf)
    
    fig = plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)
    CorrPlot([d['sumex'] for d in data],
             [d['growth_rate'] for d in data],
             [d['carbon_source'] for d in data],
             'SUMEX score',
             'Specific Growth Rate [1/h]',
             figure=fig)
    data_with_dG0 = [d for d in data if np.isfinite(d['dG0'])]
    plt.subplot(1,2,2)
    CorrPlot([-d['dG0'] for d in data_with_dG0],
             [d['growth_rate'] for d in data_with_dG0],
             [d['carbon_source'] for d in data_with_dG0],
             r'Oxidation $-\Delta_r G^\circ$ [kJ/mol]',
             'Specific Growth Rate [1/h]',
             figure=fig)
    fig.tight_layout()
    pdf.savefig(fig)

    PlotEnergy(data, 'sumex', 'SUMEX score', pdf)

    pdf.close()
    
if __name__ == "__main__":
    main()
    