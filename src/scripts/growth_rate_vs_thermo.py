import numpy as np
import csv
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.kegg import Kegg 
from pygibbs.kegg_reaction import Reaction
import matplotlib.pyplot as plt

def GetGrowthRateData():
    kegg = Kegg.getInstance()
    path = '../data/growth/growth_rates_adadi_2012.csv'
    growth_rate_dict = {}
    for row in csv.DictReader(open(path, 'r')):
        growth_rate = float(row['maximum growth rate measured'])
        carbon_source = row['carbon source']
        cid, _, _ = kegg.name2cid(carbon_source)
        if cid is None:
            raise Exception("Cannot map compound name to KEGG ID: " + carbon_source)
        growth_rate_dict[cid] = growth_rate
    return growth_rate_dict

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

def main():
    kegg = Kegg.getInstance()
    estimators = LoadAllEstimators()
    thermo = estimators['UGC']
    gr_dict = GetGrowthRateData()

    data = []
    for cid, gr in gr_dict.iteritems():
        mol = kegg.cid2mol(cid)
        mw = mol.GetExactMass()
        atom_bag, _ = mol.GetAtomBagAndCharge()
         
        r = GetFullOxidationReaction(cid)
        r.Balance(balance_water=False, balance_hydrogens=True, exception_if_unknown=True)
        dG0 = r.PredictReactionEnergy(thermo)
        s_tot = sum([abs(x) for x in r.sparse.values()])
        if np.isfinite(dG0):
            data.append((kegg.cid2name(cid), gr, dG0, atom_bag['C'], mw, s_tot))
    
    plots = []
    plots.append({'x':[-x[2] for x in data],
                  'y':[x[1] for x in data],
                  'labels':[x[0] for x in data],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/mol]',
                  'ylabel':r'Growth rate [1/hr]',
                  })

    plots.append({'x':[-x[2]/x[3] for x in data],
                  'y':[x[1] for x in data],
                  'labels':[x[0] for x in data],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/mol (C.S.)]',
                  'ylabel':r'Growth rate [1/hr]',
                  })

    plots.append({'x':[-x[2]/x[4] for x in data],
                  'y':[x[1] for x in data],
                  'labels':[x[0] for x in data],
                  'xlabel':r'Oxidation -$\Delta_r G''^\circ$ [kJ/gr]',
                  'ylabel':r'Growth rate [1/hr]',
                  })

    plots.append({'x':[x[5] for x in data],
                  'y':[x[1] for x in data],
                  'labels':[x[0] for x in data],
                  'xlabel':r'$\Sigma_i |s_i|$',
                  'ylabel':r'Growth rate [1/hr]',
                  })

    
    fig = plt.figure(figsize=(12, 12), dpi=50)

    for i, d in enumerate(plots):
        plt.subplot(2,2,i+1)
        plt.plot(d['x'], d['y'], '.', figure=fig)
        for j in xrange(len(d['x'])):
            plt.text(d['x'][j], d['y'][j], d['labels'][j], fontsize=8, figure=fig)
        plt.xlabel(d['xlabel'])
        plt.ylabel(d['ylabel'])
        r = np.corrcoef(d['x'], d['y'])
        plt.title(r'$r$ = %.2f' % r[0,1])

    
    fig.savefig('../res/growth_rate_vs_thermo.svg')
    
    
if __name__ == "__main__":
    main()