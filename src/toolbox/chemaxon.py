import logging
import subprocess
import os

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
                             calculate_all_ms=False):
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
    return diss_table, major_ms

if __name__ == "__main__":
    print "glycine", GetDissociationConstants('C(=O)(O)CN')
    print "CO2", GetDissociationConstants('O=C=O')
    print "ATP", GetDissociationConstants('Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OP(O)([O-])=O)C(O)C1O')
