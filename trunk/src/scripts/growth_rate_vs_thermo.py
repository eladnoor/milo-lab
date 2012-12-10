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
        if np.isfinite(dG0):
            data.append((kegg.cid2name(cid), gr, dG0, atom_bag['C'], mw))
    
    fig = plt.figure(figsize=(20, 7), dpi=90)
    plt.subplot(1,3,1)
    plt.plot([x[2] for x in data], [x[1] for x in data], '.', figure=fig)
    for x in data:
        plt.text(x[2], x[1], x[0], fontsize=8, figure=fig)
    plt.xlabel('Oxidation energy [kJ/mol (C.S.)]')
    plt.ylabel('Growth rate [1/hr]')

    plt.subplot(1,3,2)
    plt.plot([x[2]/x[3] for x in data], [x[1] for x in data], '.', figure=fig)
    for x in data:
        plt.text(x[2]/x[3], x[1], x[0], fontsize=8, figure=fig)
    plt.xlabel('Oxidation energy [kJ/mol of C (C.S.)]')
    plt.ylabel('Growth rate [1/hr]')

    plt.subplot(1,3,3)
    plt.plot([x[2]/x[4] for x in data], [x[1] for x in data], '.', figure=fig)
    for x in data:
        plt.text(x[2]/x[4], x[1], x[0], fontsize=8, figure=fig)
    plt.xlabel('Oxidation energy [kJ/gr (C.S.)]')
    plt.ylabel('Growth rate [1/hr]')
    
    fig.savefig('../res/growth_rate_vs_thermo.svg')
    
    
if __name__ == "__main__":
    main()