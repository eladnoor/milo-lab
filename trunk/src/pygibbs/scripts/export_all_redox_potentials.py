import sys, csv
import numpy as np
from optparse import OptionParser
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.kegg import Kegg
import collections
from pygibbs.kegg_errors import KeggParseException
from toolbox.molecule import OpenBabelError
from pygibbs.thermodynamic_constants import F, R

def CreateElementMatrix(thermo):
    kegg = Kegg.getInstance()
    atom_matrix = []
    cids = []
    for cid in thermo.get_all_cids():
        try:
            atom_vector = kegg.cid2compound(cid).get_atom_vector()
        except (KeggParseException, OpenBabelError):
            continue
        if atom_vector is not None:
            cids.append(cid)
            atom_matrix.append(atom_vector)
    atom_matrix = np.array(atom_matrix)
    return cids, atom_matrix

def FindRedoxPairs(cids, atom_matrix, thermo):
    """
        Finds pairs of rows in the matrix where the only difference
        is exactly one or two electrons (H is ignored)
    """
    def hash_atom_vector(v):
        return ','.join(["%d:%d" % (i, v[i])
                         for i in np.nonzero(v)[0]
                         if i != 1]) # ignore H atoms
    
    dG0_f = thermo.GetTransformedFormationEnergies(cids)

    kegg = Kegg.getInstance()

    # remove the H column, and create two matrices with one or two added
    # electrons
    lookup = collections.defaultdict(list)
    for i in xrange(atom_matrix.shape[0]):
        h = hash_atom_vector(atom_matrix[i, :])
        lookup[h].append(i)
    
    rowdicts = []
    fieldnames = ['name_ox', 'CID_ox', 'ne_ox', 'name_red', 'CID_red', 'ne_red',
                  'E_tag', 'pH', 'ref']
    for delta_e in [1, 2]:
        # check if there are rows which equal other rows, with a difference
        # of 'delta_e' electrons
        reduction_vector = np.zeros((atom_matrix.shape[1]), dtype='int')
        reduction_vector[0] = delta_e
        for i in xrange(atom_matrix.shape[0]):
            h = hash_atom_vector(atom_matrix[i, :] + reduction_vector)
            partners = lookup[h]
            for k in partners:
                #delta_H = atom_matrix[k, 1] - atom_matrix[i, 1]
                ddG0 = dG0_f[i, 0] - dG0_f[k, 0]
                E_prime = ddG0 / (F * delta_e)
                rowdicts.append({
                    "name_ox":kegg.cid2name(cids[i]),
                    "CID_ox":cids[i],
                    "ne_ox":atom_matrix[i, 0],
                    "name_red":kegg.cid2name(cids[k]),
                    "CID_red":cids[k],
                    "ne_red":atom_matrix[k, 0],
                    "E_tag":"%.2f" % E_prime,
                    "pH":"%g" % thermo.pH,
                    "ref":thermo.name})
    return rowdicts, fieldnames

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="PRC",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-c", "--csv_out_fname",
                          dest="csv_out_fname",
                          default="../res/redox_pair.csv",
                          help="CSV output filename")
    return opt_parser

def ExportThermo():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    thermo = estimators[options.thermodynamics_source]

    print 'Thermodynamic source:', thermo.name
    print 'CSV output filename:', options.csv_out_fname

    cids, atom_matrix = CreateElementMatrix(thermo)
    rowdicts, fieldnames = FindRedoxPairs(cids, atom_matrix, thermo)
    csv_out = csv.DictWriter(open(options.csv_out_fname, 'w'), fieldnames)
    csv_out.writeheader()
    csv_out.writerows(rowdicts)
    
if __name__ == "__main__":
    ExportThermo()
