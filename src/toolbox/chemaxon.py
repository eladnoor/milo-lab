import logging, itertools, subprocess, os
import numpy as np

CXCALC_BIN = "/usr/bin/cxcalc"

class ChemAxonError(Exception):
    pass

def RunCxcalc(args):
    if not os.path.exists(CXCALC_BIN):
        raise Exception("Jchem must be installed to calculate pKa data.")
    
    debug_args = "ARGS: %s" % ' '.join([CXCALC_BIN] + args)
    logging.debug("\n" + debug_args)
    p = subprocess.Popen([CXCALC_BIN] + args,
                         executable=CXCALC_BIN, stdout=subprocess.PIPE)
    #p.wait()
    #os.remove(temp_fname)
    res = p.communicate()[0]
    if p.returncode != 0:
        raise ChemAxonError(debug_args)
    logging.debug("OUTPUT: %s" % res)
    return res

def ParsePkaOutput(s, n_acidic, n_basic, pH_list):
    """
        Returns:
            A dictionary that maps the atom index to a list of pKas
            that are assigned to that atom.
    """
    atom2pKa = {}

    pkaline = s.split('\n')[1]
    splitline = pkaline.split('\t')
    splitline.pop(0)
    
    if n_acidic + n_basic > 0:
        if len(splitline) != (n_acidic + n_basic + 1 + len(pH_list)):
            raise ChemAxonError('ChemAxon failed to find any pKas')
        
        pKa_list = []
        acid_or_base_list = []
        for i in range(n_acidic + n_basic):
            x = splitline.pop(0) 
            if x == '':
                continue
            
            pKa_list.append(float(x))
            if i < n_acidic:
                acid_or_base_list.append('acid')
            else:
                acid_or_base_list.append('base')
        
        atom_list = splitline.pop(0)

        if atom_list: # a comma separated list of the deprotonated atoms
            atom_numbers = [int(x)-1 for x in atom_list.split(',')]
            for i, j in enumerate(atom_numbers):
                atom2pKa.setdefault(j, [])
                atom2pKa[j].append((pKa_list[i], acid_or_base_list[i]))
    
    smiles_list = splitline
    return atom2pKa, smiles_list

def _GetDissociationConstants(molstring, n_acidic=10, n_basic=10, pH_list=None):
    """
        Returns:
            A pair of (pKa list, major pseudoisomer)
            
            - the pKa list is of the pKa values in ascending order.
            - the major pseudoisomer is a SMILES string of the major species
              at the given pH.
    """
    pH_list = pH_list or []
    args = []
    if n_acidic + n_basic > 0:
        args += ['pka', '-a', str(n_acidic), '-b', str(n_basic),
                 '-M', 'true', '-P', 'dynamic']
    for pH in pH_list:
        args += ['majorms', '--pH', str(pH), '-M', 'true']
    args += [molstring]
    
    output = RunCxcalc(args)
    atom2pKa, smiles_list = ParsePkaOutput(output, n_acidic, n_basic, pH_list)
    
    all_pKas = []
    for pKa_list in atom2pKa.values():
        all_pKas += [pKa for pKa, _ in pKa_list]
    
    return sorted(all_pKas), smiles_list

def GetDissociationConstants(molstring, n_acidic=10, n_basic=10, mid_pH=7,
                             calculate_all_ms=False, transform_multiples=True):
    """
        Arguments:
            molstring - a text description of the molecule (SMILES or InChI)
            n_acidic  - the max no. of acidic pKas to calculate
            n_basic   - the max no. of basic pKas to calculate
            mid_pH    - the default pH for which the major pseudoisomer is calculated
            calculate_all_ms - if True, uses ChamAxon to find all SMILES descriptors
                               of all pseudoisomers. Otherwise only the major one
                               is returned.
        
        Returns a pair:
            (diss_constants, major_ms)
            
        - diss_constants is a list of 3-tuples: (pKa, SMILES below, SMILES above)
        - major_ms is a SMILES string of the major pseudoisomer at pH_mid 
    """
    all_pKas, smiles_list = _GetDissociationConstants(molstring, n_acidic, 
                                                    n_basic, [mid_pH])
    major_ms = smiles_list[0]
    if all_pKas == []:
        return [], major_ms
    
    if calculate_all_ms:
        pH_list = [all_pKas[0]-0.1]
        for i in range(len(all_pKas)-1):
            pH_list.append((all_pKas[i] + all_pKas[i+1]) * 0.5)
        pH_list.append(all_pKas[-1]+0.1)
        
        _, smiles_list = _GetDissociationConstants(molstring, 0, 0, pH_list)
    else: # find the index of the major pseudoisomer and use None for all the others
        smiles_list = [None] * (len(all_pKas) + 1)
        if mid_pH < min(all_pKas):
            major_i = 0
        else:
            major_i = 1 + max([i for i, pKa in enumerate(all_pKas) if pKa < mid_pH])
        smiles_list[major_i] = major_ms
    
    diss_table = []
    for i, pKa in enumerate(all_pKas):
        diss_table.append((pKa, smiles_list[i], smiles_list[i+1]))
    
    if transform_multiples:
        diss_table= _TransformMultiples(diss_table)
    return diss_table, major_ms

def _TransformMultiples(diss_table):
    """
        There are two ways to interpret the pKa values coming from ChemAxon.
        One way is to say each one corresponds to a specific protonation site,
        and represents the equilibrium constant for that site (regardless of the
        environment).
        The other way is to say each pKa is the equilibrium between the family 
        of species with some nH and the family of species with nH+1.
        
        This method will take ChemAxon results, and recalculate the pKas assuming
        the first interpretation is correct. The return value will look exactly 
        the same as the input (i.e. a diss_table) but the pKas will
        now represent the sum of all pseudoisomers in the relevant charge groups.
    """
    # sort in descending order of pKa, i.e. ascending dG0 (which corresponds to
    # ascending nH). The first value in this vector will be the difference in 
    # dG0 between the species with the smallest nH and nH+1. 
    diss_table = sorted(diss_table, key=lambda x: -x[0])
    
    # Ka[i] = [i][H+]/[i+1]
    Ka_list = [10**(-pKa) for (pKa, _, _) in diss_table]

    relative_conc = [1]
    # In this case, we must calculate the ddG0 for all pseudoisomers with the
    # same nH (i.e. combinatorically using all protonations states which sum
    # up to a specific value). The product of the Ka values is the ratio between
    # that species and the fully deprotonated one.
    for i in xrange(len(Ka_list)):
        sum_conc = 0
        for Ka_subset in itertools.combinations(Ka_list, i+1): # all choices of i values from the Ka list
            sum_conc += np.prod(Ka_subset)
        relative_conc.append(sum_conc)
        Ka_i = relative_conc[i+1] / relative_conc[i]
        (_, smiles1, smiles2) = diss_table[i]
        diss_table[i] = (-np.log10(Ka_i), smiles1, smiles2)
    
    return diss_table

if __name__ == "__main__":
    
    diss_table_example = [(4.0,None,None), (4.01, None, None), (9, None, None)]
    new_diss_table = _TransformMultiples(diss_table_example)
    
    print diss_table_example
    print new_diss_table
    
    from toolbox.molecule import Molecule
    compound_list = [('glycine', 'C(=O)(O)CN'),
                     ('CO2', 'O=C=O'),
                     ('ATP', 'Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OP(O)([O-])=O)C(O)C1O'),
                     ('3-Ketoarabinitol', 'OCC(O)C(C(O)CO)=O')]
    
    for name, smiles in compound_list:
        diss_table, major_ms = GetDissociationConstants(smiles, transform_multiples=False)
        diss_table2 = _TransformMultiples(diss_table)
        m = Molecule.FromSmiles(major_ms)
        print name, m.ToInChI()
        for i in xrange(len(diss_table)):
            print diss_table[i][0], diss_table2[i][0]
