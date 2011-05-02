#!/usr/bin/python

import csv
import logging
from toolbox.molecule import Molecule

class PseudoisomerEntry(object):
    def __init__(self, net_charge, hydrogens, magnesiums, smiles="",
                 dG0=None, cid=None, name=None, ref='', use_for=None):
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
    
    def Clone(self):
        return PseudoisomerEntry(self.net_charge, self.hydrogens, self.magnesiums,
            self.smiles, self.dG0, self.cid, self.name, self.ref, self.use_for)
    
    def Complete(self):
        """Returns true if it has enough data to use for training/testing."""
        if not self.name or not self.smiles:
            return False
        
        try:
            int(self.net_charge)
            int(self.hydrogens)
            int(self.magnesiums)
        except Exception:
            return False
        
        return True
        
    def Train(self):
        return self.use_for.lower() == 'train'
    
    def Test(self):
        return self.use_for.lower() == 'test'
    
    def Skip(self):
        return self.use_for.lower() == 'skip'
    
    def Mol(self):
        """Returns a new Mol object corresponding to this compound."""
        if not self.smiles:
            return None
        return Molecule.FromSmiles(self.smiles)
    
    def MolNoH(self):
        mol = self.Mol()
        if not mol:
            return None
        
        mol.RemoveHydrogens()
        return mol
    
    def __str__(self):
        return '%30s (z=%2d, nH=%2d, nMg=%2d): dG0=%7.1f' % \
            (self.name, self.net_charge or 0, self.hydrogens or 0, 
             self.magnesiums or 0, self.dG0)
        
    def __hash__(self):
        return hash((self.name, self.net_charge,
                     self.hydrogens, self.magnesiums))
    
    def Tag(self):
        return '%s%d' % (self.name, hash(self))
        
    tag = property(Tag)
    
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
            name = row_dict.get('name')
            ref = row_dict.get('ref')
            smiles = row_dict.get('smiles')
            use_for = row_dict.get('use for')
            charge = row_dict.get('z')
            hydrogens = row_dict.get('nH')
            nMg = row_dict.get('nMg')
            cid = row_dict.get('cid')
            dG0 = row_dict.get('dG0')
            
            if not charge:
                logging.warning('Failed to read charge for compound %s', name)
                charge = None
            else:
                charge = int(charge)
            
            if not hydrogens:
                logging.warning('Failed to read hydrogens for compound %s', name)
                hydrogens = None
            else:
                hydrogens = int(hydrogens)
            
            if not nMg:
                logging.warning('Failed to read magnesiums for compound %s', name)
                nMg = None
            else:
                nMg = int(nMg)
            
            if not cid:
                logging.warning('Failed to read KEGG ID for compound %s', name)
                cid = None
            else:
                cid = int(cid)
            
            if not dG0:
                logging.warning('Failed to read dG0 for compound %s', name)
                dG0 = None
            else:
                dG0 = float(dG0)
                
            if use_for == 'skip':
                continue
            
            comp = PseudoisomerEntry(charge, hydrogens, nMg, smiles,
                                     dG0=dG0, name=name, cid=cid,
                                     ref=ref, use_for=use_for)
                
            logging.info('Reading data for %s (C%05d)' % (comp, cid or -1))
            
            pseudoisomers.append(comp)
        
        return PseudoisomersData(pseudoisomers)

    def ReadDissociationData(self, dissociation):
        new_pseudoisomers = {}
        for pdata in self.pseudoisomers:
            cid = pdata.cid
            nH = pdata.hydrogens
            nMg = pdata.magnesiums
            ref = pdata.ref
            if not cid or nH == None: # it is important not to use ("or not nH") since it can be a 0
                continue
            if nMg == None:
                nMg = 0
            new_pseudoisomers[cid, nMg, nH, ref] = pdata
            pK_table = dissociation.GetDissociationTable(cid)
            if pK_table == None:
                continue
            try:
                for new_pdata in pK_table.GenerateAllPseudoisomerEntries(pdata):
                    cid = new_pdata.cid
                    nH = new_pdata.hydrogens
                    nMg = new_pdata.magnesiums
                    ref = new_pdata.ref
                    if (cid, nMg, nH) not in new_pseudoisomers:
                        new_pseudoisomers[cid, nMg, nH, ref] = new_pdata
            except KeyError, msg:
                logging.error(msg)
                raise Exception('Cannot find the pK data for compound C%05d' % cid)

        species_set = set()                
        for cid, nMg, nH, ref, pdata in sorted(new_pseudoisomers.iteritems()):
            if (cid, nMg, nH) in species_set:
                raise Exception("Conflict: " + str([cid, nMg, nH]))
            else:
                species_set.add((cid, nMg, nH))
            print pdata, ref

if __name__ == '__main__':
    from pygibbs.dissociation_constants import DissociationConstants

    pdata = PseudoisomersData.FromFile(
        '../data/thermodynamics/dG0.csv')
    dissociation = DissociationConstants.FromFile(
        '../data/thermodynamics/dissociation_constants.csv')
    pdata.ReadDissociationData(dissociation)
    