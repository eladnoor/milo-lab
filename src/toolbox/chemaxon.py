import logging
import subprocess
import os
from pygibbs.thermodynamic_constants import default_pH

CXCALC_BIN = "/home/eladn/opt/jchem-5.5.1.0/bin/cxcalc"

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
    pH2smiles = {}

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
    
    for pH in pH_list:
        pH2smiles[pH] = splitline.pop(0)
    
    return atom2pKa, pH2smiles

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
    atom2pKa, pH2smiles = ParsePkaOutput(output, n_acidic, n_basic, pH_list)
    
    all_pKas = []
    for pKa_list in atom2pKa.values():
        all_pKas += [pKa for pKa, _ in pKa_list]
    
    return sorted(all_pKas), pH2smiles

def GetDissociationConstants(molstring, n_acidic=10, n_basic=10):
    all_pKas, pH2smiles = _GetDissociationConstants(molstring, n_acidic, 
                                                    n_basic, [])
    if all_pKas == []:
        return all_pKas, pH2smiles
    
    pH_list = [all_pKas[0]-0.1]
    for i in range(len(all_pKas)-1):
        pH_list.append((all_pKas[i] + all_pKas[i+1]) * 0.5)
    pH_list.append(all_pKas[-1]+0.1)
    
    _, pH2smiles = _GetDissociationConstants(molstring, 0, 0, pH_list)
    
    return all_pKas, pH2smiles

if __name__ == "__main__":
    pH_list = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0]
    print "urea", GetDissociationConstants('C(=O)(N)N')
    print "glycine", GetDissociationConstants('C(=O)(O)CN')
    print "ATP", GetDissociationConstants('Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O')