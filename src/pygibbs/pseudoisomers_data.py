#!/usr/bin/python

import csv
import logging
import pybel
from pygibbs.thermodynamic_constants import R, default_T, dG0_f_Mg, debye_huckel
import pylab
from toolbox.util import log_sum_exp
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.kegg import Kegg

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
        except Exception, e:
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
        return pybel.readstring('smiles', self.smiles)
    
    def MolNoH(self):
        mol = self.Mol()
        if not mol:
            return None
        
        mol.removeh()
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
    
class DissociationTable(object):
    
    def __init__(self, cid=None):
        self.ddGs = {}
        self.cid = cid
        self.min_nH = None # the nH of the most basic pseudoisomer
        self.min_charge = None # the charge of the most basic pseudoisomer
        self.min_dG0 = 0 # the dG0 of the most basic pseudoisomer
        self.kegg = Kegg.getInstance()

    @staticmethod
    def ReadDissociationCsv(filename='../data/thermodynamics/dissociation_constants.csv'):
        """
            Parses a CSV file that contains pKa and pKMg data for many compounds
            and returns a dictionary of their DissociationTables, where the key
            is the CID.
        """
        cid2pK = {}

        csv_reader = csv.DictReader(open(filename, 'r'))
        for i, row in enumerate(csv_reader):
            if not row['cid']:
                continue # without a CID we cannot match this to the dG0 table
            cid = int(row['cid'])
            logging.debug("Parsing row #%d, compound C%05d" % (i, cid))

            nH_below = int(row['nH_below'])
            nH_above = int(row['nH_above'])
            nMg_below = int(row['nMg_below'])
            nMg_above = int(row['nMg_above'])
            ref = row['ref']
            T = float(row['T'] or default_T)
            cid2pK.setdefault(cid, DissociationTable(cid))
            cid2pK[cid].min_nH = min(nH_above, cid2pK[cid].min_nH or nH_above)


            if row['type'] == 'acid-base':
                pKa = float(row['pK'])
                if nMg_below != nMg_above:
                    raise Exception('C%05d has different nMg below and above '
                                    'the pKa = %.1f' % pKa)
                cid2pK[cid].AddpKa(pKa, nH_below, nH_above, nMg_below, ref, T)

            if row['type'] == 'Mg':
                pKMg = float(row['pK'])
                if nH_below != nH_above:
                    raise Exception('C%05d has different nH below and above '
                                    'the pK_Mg = %.1f' % pKMg)
                cid2pK[cid].AddpKMg(pKMg, nMg_below, nMg_above, nH_below, ref, T)
        
        for pK_table in cid2pK.values():
            pK_table.CalculateCharge()
        
        return cid2pK
    
    def AddpKa(self, pKa, nH_below, nH_above, nMg=0, ref="", T=default_T):
        if nH_below != nH_above+1:
            raise Exception('A H+ dissociation constant (pKa) has to represent an '
                            'increase of exactly one hydrogen')
        
        ddG0 = R * T * pylab.log(10) * pKa
        self.ddGs[(nH_above, nH_below, nMg, nMg)] = (-ddG0, ref) # adding H+ decreases dG0
        
    def AddpKMg(self, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T):
        if nMg_below != nMg_above+1:
            raise Exception('A Mg+2 dissociation constant (pK_Mg) has to represent an '
                            'increase of exactly one magnesium ion')

        ddG0 = R * T * pylab.log(10) * pKMg - dG0_f_Mg
        self.ddGs[(nH, nH, nMg_above, nMg_below)] = (-ddG0, ref) # adding Mg+2 decreases dG0
    
    def GetSingleStep(self, nH_from, nH_to, nMg_from, nMg_to):
        if nH_from == nH_to and nMg_from == nMg_to:
            return 0, None
        
        if nH_from != nH_to and nMg_from != nMg_to:
            raise Exception('A dissociation constant can either represent a'
                ' change in hydrogens or in magnesiums, but not both')
        
        try:
            if nMg_to == nMg_from+1 or nH_to == nH_from+1:
                ddG0, ref = self.ddGs[nH_from, nH_to, nMg_from, nMg_to]
                return (ddG0, ref)
            if nMg_from == nMg_to+1:
                ddG0, ref = self.ddGs[nH_from, nH_to, nMg_to, nMg_from]
                return (-ddG0, ref)
            if nH_from == nH_to+1:
                ddG0, ref = self.ddGs[nH_to, nH_from, nMg_from, nMg_to]
                return (-ddG0, ref)
        except KeyError:
            raise KeyError('The dissociation constant for (nH=%d,nMg=%d) -> '
                            '(nH=%d,nMg=%d) is missing' % (nH_from, nMg_from, nH_to, nMg_to))

        raise Exception('A dissociation constant can either represent a'
                ' change in only one hydrogen or magnesium')
    
    def ConvertPseudoisomer(self, dG0, nH_from, nH_to=None, nMg_from=0, nMg_to=0):
        if not nH_to:
            nH_to = self.min_nH
        
        pdata = PseudoisomerEntry(net_charge=0, hydrogens=nH_from, 
            magnesiums=nMg_from, smiles="", dG0=dG0)
        comp = self.ConvertPseudoisomerEntry(pdata, nH_to, nMg_to)
        return comp.dG0
    
    def ConvertPseudoisomerEntry(self, pdata, nH_to, nMg_to):
        """
            Returns the difference in dG0 between any two pseudoisomers.
        """
        nH_from = pdata.hydrogens
        nMg_from = pdata.magnesiums
        
        step_list = []
        
        # first remove all Mgs from the original
        step_list += [self.GetSingleStep(nH_from, nH_from, nMg, nMg-1)
                      for nMg in xrange(nMg_from, 0, -1)]

        # then change the nH to fit the target (using only species without Mg)
        if nH_from < nH_to:
            step_list += [self.GetSingleStep(nH, nH+1, 0, 0)
                          for nH in xrange(nH_from, nH_to)]
        elif nH_from > nH_to:
            step_list += [self.GetSingleStep(nH, nH-1, 0, 0)
                          for nH in xrange(nH_from, nH_to, -1)]
        
        # finally add back all the Mg in the target
        step_list += [self.GetSingleStep(nH_to, nH_to, nMg, nMg+1)
                      for nMg in xrange(0, nMg_to)]
        
        
        total_ddG0 = sum([ddG0 for ddG0, _ref in step_list])
        total_ref = ';'.join([(ref or "") for _ddG0, ref in step_list])
        
        comp = pdata.Clone()
        comp.dG0 += total_ddG0
        comp.ref += ';' + total_ref
        comp.smiles = ''
        comp.hydrogens = nH_to
        comp.magnesiums = nMg_to
        comp.net_charge += (nH_to - nH_from) + 2 * (nMg_to - nMg_from)
        
        return comp
    
    def CalculateCharge(self):
        # get the charge and nH of the default pseudoisomer in KEGG:
        z = self.kegg.cid2charge(self.cid, correctForPH=False) or 0
        nH = self.kegg.cid2num_hydrogens(self.cid, correctForPH=False) or 0
        
        # calculate the charge for the most basic species
        self.min_charge = z + (self.min_nH - nH)
        
    def GenerateAll(self):
        if self.min_charge == None:
            raise Exception('The minimal charge has to be set before generating'
            ' all psuedoisomers')
        if self.min_dG0 == None:
            raise Exception('The base formation energy has to be set before generating'
            ' all psuedoisomers')
                
        pdata = PseudoisomerEntry(net_charge=self.min_charge, hydrogens=self.min_nH,
            magnesiums=0, smiles="", dG0=self.min_dG0)
        return self.GenerateAllPseudoisomerEntries(pdata)
    
    def GenerateAllPseudoisomerEntries(self, pdata):
        pseudoisomers = {}
        pseudoisomers[pdata.hydrogens, pdata.magnesiums] = pdata.Clone()

        for nH_above, nH_below, nMg_above, nMg_below in self.ddGs.keys():
            
            if (nH_above, nMg_above) not in pseudoisomers:
                pseudoisomers[nH_above, nMg_above] = \
                    self.ConvertPseudoisomerEntry(pdata, nH_above, nMg_above)

            if (nH_below, nMg_below) not in pseudoisomers:
                pseudoisomers[nH_below, nMg_below] = \
                    self.ConvertPseudoisomerEntry(pdata, nH_below, nMg_below)
        
        return pseudoisomers.values()
    
    def Transform(self, pH, I, pMg, T):
        # assume that the dG0_f of the most basic psuedoisomer is 0, 
        # and calculate the transformed dG'0_f relative to it.
        
        dG0_tag_vec = []
        for pseudoisomer in self.GenerateAll():
            nH = pseudoisomer.hydrogens
            nMg = pseudoisomer.magnesiums
            z = pseudoisomer.net_charge
            dG0 = pseudoisomer.dG0
            
            DH = debye_huckel(I)
            dG0_tag = dG0 + \
                      nMg * (R*T*pylab.log(10)*pMg - dG0_f_Mg) + \
                      nH  * (R*T*pylab.log(10)*pH + DH) - (z**2) * DH
            dG0_tag_vec.append(dG0_tag)
        
        dG0_tag_total = -R * T * log_sum_exp([g / (-R*T) for g in dG0_tag_vec])
        
        return dG0_tag_total
        
    def GetPseudoisomerMap(self):
        pmap = PseudoisomerMap()
        for pdata in self.GenerateAll():
            pmap.Add(nH=pdata.hydrogens, z=pdata.net_charge, 
                     nMg=pdata.magnesiums, dG0=pdata.dG0)
        return pmap

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

    def ReadDissociationData(self):
        cid2pK = DissociationTable.ReadDissociationCsv()
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
            if cid not in cid2pK:
                continue
            
            pK_table = cid2pK[cid]
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
        for (cid, nMg, nH, ref), pdata in sorted(new_pseudoisomers.iteritems()):
            if (cid, nMg, nH) in species_set:
                raise Exception("Conflict: " + str([cid, nMg, nH]))
            else:
                species_set.add((cid, nMg, nH))
            print pdata, ref

if __name__ == '__main__':
    pdata = PseudoisomersData.FromFile('../data/thermodynamics/dG0.csv')
    pdata.ReadDissociationData()
    