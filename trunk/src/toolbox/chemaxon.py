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
    
    logging.debug("\nARGS: %s" % ' '.join([CXCALC_BIN] + args))
    p = subprocess.Popen([CXCALC_BIN] + args,
                         executable=CXCALC_BIN, stdout=subprocess.PIPE)
    #p.wait()
    #os.remove(temp_fname)
    res = p.communicate()[0]
    if p.returncode != 0:
        raise ChemAxonError()
    logging.debug("OUTPUT: %s" % res)
    return res

def ParsePkaOutput(s, n_acidic, n_basic):
    """
        Returns:
            A dictionary that maps the atom index to a list of pKas
            that are assigned to that atom.
    """
    pkaline = s.split('\n')[1]
    splitline = pkaline.split('\t')
    if len(splitline) != (n_acidic + n_basic + 3):
        raise ChemAxonError('ChemAxon failed to find any pKas')
    
    apKa_list = [float(x) for x in splitline[1:(n_acidic+1)] if x != '']
    bpKa_list = [float(x) for x in splitline[(n_acidic+1):(n_acidic+n_basic+1)] if x != '']
             
    pKa_list = apKa_list + bpKa_list
    acid_or_base_list = ['acid'] * len(apKa_list) + ['base'] * len(bpKa_list)        
       
    [atom_list, smiles] = splitline[n_acidic+n_basic+1:n_acidic+n_basic+3] 

    atom2pKa = {}
    if atom_list: # a comma separated list of the deprotonated atoms
        atom_numbers = [int(x)-1 for x in atom_list.split(',')]
        for i, j in enumerate(atom_numbers):
            atom2pKa.setdefault(j, [])
            atom2pKa[j].append((pKa_list[i], acid_or_base_list[i]))
    
    if smiles: # a SMILES string of the major microspecies
        return atom2pKa, smiles
    else:
        return atom2pKa, None

def GetDissociationConstants(molstring, n_acidic=10, n_basic=10, pH=default_pH):
    """
        Returns:
            A pair of (pKa list, major pseudoisomer)
            
            - the pKa list is a list of the pKa values in ascending order.
            - the major pseudoisomer is a SMILES string of the major species
              at the given pH.
    """
    args = ['pka',
            '-a', str(n_acidic),
            '-b', str(n_basic),
            '-M', 'true',
            '-P', 'dynamic',
            'majorms',
            '--pH', str(pH),
            molstring]
    
    output = RunCxcalc(args)
    atom2pKa, smiles = ParsePkaOutput(output, n_acidic, n_basic)
    
    all_pKas = []
    for pKa_list in atom2pKa.values():
        all_pKas += [pKa for pKa, _ in pKa_list]
    
    return sorted(all_pKas), smiles

if __name__ == "__main__":
    print "urea", GetDissociationConstants('C(=O)(N)N')
    print "glycine", GetDissociationConstants('C(=O)(O)CN')
    print "ATP", GetDissociationConstants('Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O')