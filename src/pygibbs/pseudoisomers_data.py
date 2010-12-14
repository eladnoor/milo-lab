#!/usr/bin/python

import csv
import logging
import pybel


class PseudoisomerEntry(object):
    def __init__(self, net_charge, hydrogens, magnesiums, smiles,
                 dG0=None, cid=None, name=None, ref=None, use_for=None):
        """Initialize a compound."""
        self.cid = cid
        self.name = name
        self.net_charge = net_charge
        self.hydrogens = hydrogens
        self.magnesiums = magnesiums
        self.smiles = smiles
        self.dG0 = dG0
        self.ref = ref
        self.use_for = use_for
        
    def Train(self):
        return self.use_for.lower() == 'train'
    
    def Test(self):
        return self.use_for.lower() == 'test'
    
    def Skip(self):
        return self.use_for.lower() == 'skip'
    
    def Mol(self):
        """Returns a new Mol object corresponding to this compound."""
        return pybel.readstring('smiles', self.smiles)
    
    def __str__(self):
        return '%s (z=%s, nH=%s, nMg=%s)' % (self.name, self.net_charge,
                                             self.hydrogens, self.magnesiums)
    

class PseudoisomersData(object):
    
    def __init__(self, pseudoisomers):
        """Initialize PseudoisomersData."""
        self.pseudoisomers = pseudoisomers
    
    def __iter__(self):
        return iter(self.pseudoisomers)
    
    @staticmethod
    def FromFile(filename):
        """Build a PseudoisomersData object from a file."""
        pseudoisomers = []
        for row_dict in csv.DictReader(open(filename)):
            name = row_dict.get('compound name')
            ref = row_dict.get('ref')
            smiles = row_dict.get('smiles')
            use_for = row_dict.get('use for')
            charge = row_dict.get('charge')
            hydrogens = row_dict.get('hydrogens')
            nMg = row_dict.get('Mg')
            cid = row_dict.get('cid')
            dG0 = row_dict.get('dG0')
            
            if not charge:
                logging.warning('Failed to read charge for compound %s', name)
            else:
                charge = int(charge)
            
            if not hydrogens:
                logging.warning('Failed to read hydrogens for compound %s', name)
            else:
                hydrogens = int(hydrogens)
            
            if not nMg:
                logging.warning('Failed to read magnesiums for compound %s', name)
            else:
                nMg = int(nMg)
            
            if not cid:
                logging.warning('Failed to read KEGG ID for compound %s', name)
            else:
                cid = int(cid)
            
            if not dG0:
                logging.warning('Failed to read dG0 for compound %s', name)
            else:
                dG0 = float(dG0)
            
            comp = PseudoisomerEntry(charge, hydrogens, nMg, smiles,
                                     dG0=dG0, name=name, cid=cid,
                                     ref=ref, use_for=use_for)
            logging.info('Reading data for %s', comp)
            
            pseudoisomers.append(comp)
        
        return PseudoisomersData(pseudoisomers)


if __name__ == '__main__':
    pdata = PseudoisomersData.FromFile('../data/thermodynamics/dG0.csv')
    for pisomer in pdata:
        print pisomer
    