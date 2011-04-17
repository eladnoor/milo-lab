import csv, logging
from kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamic_constants import R, default_T, dG0_f_Mg, debye_huckel
import numpy as np
from toolbox.util import log_sum_exp
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.pseudoisomers_data import PseudoisomerEntry

class DissociationConstants(object):
    
    def __init__(self):
        self.cid2DissociationTable = {}
    
    @staticmethod
    def FromFile(filename='../data/thermodynamics/dissociation_constants.csv'):
        """
            Parses a CSV file that contains pKa and pKMg data for many compounds
            and returns a dictionary of their DissociationTables, where the key
            is the CID.
        """
        diss = DissociationConstants()

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
            smiles_below = row['smiles_below']
            smiles_above = row['smiles_above']
            
            ref = row['ref']
            T = float(row['T'] or default_T)
            diss.SetMinNumHydrogens(cid, nH_above)

            if row['type'] == 'acid-base':
                pKa = float(row['pK'])
                if nMg_below != nMg_above:
                    raise Exception('C%05d has different nMg below and above '
                                    'the pKa = %.1f' % (cid, pKa))
                diss.AddpKa(cid, pKa, nH_below, nH_above, nMg_below, ref, T, smiles_below, smiles_above)
            elif row['type'] == 'Mg':
                pKMg = float(row['pK'])
                if nH_below != nH_above:
                    raise Exception('C%05d has different nH below and above '
                                    'the pK_Mg = %.1f' % pKMg)
                try:
                    diss.AddpKMg(cid, pKMg, nMg_below, nMg_above, nH_below, ref, T, smiles_below, smiles_above)
                except Exception, e:
                    raise Exception("In C%05d: %s" % (cid, str(e)))
            elif row['pK']:
                raise ValueError('The row about C%05d has a pK although it is not "acid-base" nor "Mg"' % cid)
            elif nMg_below != nMg_above:
                raise ValueError('The row about C%05d has different nMgs although it is not "Mg"' % cid)
            elif nH_below != nH_above:
                raise ValueError('The row about C%05d has different nHs although it is not "acid-base"' % cid)
        
        diss.CalculateAllCharges()
        return diss
    
    def AddpKa(self, cid, pKa, nH_below, nH_above, nMg=0, ref="", T=default_T, 
               smiles_below=None, smiles_above=None):
        diss_table = self.GetDissociationTable(cid)
        diss_table.AddpKa(pKa, nH_below, nH_above, nMg,
                          ref, T, smiles_below, smiles_above)
        
    def AddpKMg(self, cid, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T, 
                smiles_below=None, smiles_above=None):
        diss_table = self.GetDissociationTable(cid)
        diss_table.AddpKMg(pKMg, nMg_below, nMg_above, nH, 
                          ref, T, smiles_below, smiles_above)
        
    def SetMinNumHydrogens(self, cid, min_nH):
        """
            Sets the minimal number of hydrogen atoms for a specified CID.
            
            If the value is not lower than the previously provided min_nH, the
            lower value is kept unchanged. 
        """
        diss_table = self.GetDissociationTable(cid)
        if not diss_table.min_nH or diss_table.min_nH > min_nH: 
            diss_table.min_nH = min_nH
        
    def GetDissociationTable(self, cid):
        if cid not in self.cid2DissociationTable:
            diss = DissociationTable(cid)
            self.cid2DissociationTable[cid] = diss
            return diss
        else:
            return self.cid2DissociationTable[cid]
    
    def CalculateAllCharges(self):
        for diss_table in self.cid2DissociationTable.values():
            diss_table.CalculateCharge()
    
    def ToDatabase(self, db, table_name):
        """
            Load the data regarding pKa values according to KEGG compound IDs.
        """
        
        db.CreateTable(table_name, """
            cid INT, name TEXT, 
            nH_below INT, nH_above INT, 
            nMg_below INT, nMg_above INT, 
            smiles_below TEXT, smiles_above TEXT, 
            ddG REAL, ref TEXT""")
        for diss_table in self.cid2DissociationTable.values():
            diss_table.ToDatabase(db, table_name)

        db.Commit()

    def GetAllCids(self):
        return set(self.cid2DissociationTable.keys())

    def ReverseTranformNistRows(self, nist_rows):
        kegg = Kegg.getInstance()

        encountered_cids = self.GetAllCids()

        # For each CID which doesn't have a known pKa or pKMg, assume that
        # that there are none in the relevant range of pH and pMg.
        # Therefore, add an empty DissociationTable to each one of them.
        for nist_row_data in nist_rows:
            cids_in_reaction = nist_row_data.GetAllCids()
            for cid in cids_in_reaction.difference(encountered_cids):
                encountered_cids.add(cid)
                try:
                    min_nH, min_charge = kegg.cid2nH_and_charge(cid)
                    diss = self.GetDissociationTable(cid) # this creates an empty table
                    diss.min_nH, diss.min_charge = min_nH, min_charge
                except TypeError:
                    logging.warning('cannot add %s (C%05d) since nH or charge '
                                    'cannot be determined' % 
                                    (kegg.cid2name(cid), cid))
            
        all_cids_with_pKa = self.GetAllCids()

        data = {}
        data['cids_to_estimate'] = sorted(all_cids_with_pKa)
        
        # the transformed (observed) free energy of the reactions dG'0_r
        data['dG0_r_tag'] = np.zeros((0, 1))
        
        # dG'0_r - dG0_r  (which is only a function of the conditions and pKas)
        data['ddG0_r'] = np.zeros((0, 1))
        
        data['pH'] = np.zeros((0, 1))
        data['I'] = np.zeros((0, 1))
        data['pMg'] = np.zeros((0, 1))
        data['T'] = np.zeros((0, 1))
        data['S'] = np.zeros((0, len(data['cids_to_estimate']))) # stoichiometric matrix
        
        for nist_row_data in nist_rows:
            # check that all participating compounds have a known pKa
            cids_in_reaction = nist_row_data.GetAllCids()
            cids_without_pKa = cids_in_reaction.difference(all_cids_with_pKa)
            if cids_without_pKa:
                logging.debug('reaction contains CIDs with unknown pKa values: %s' % \
                              ', '.join(['C%05d' % cid for cid in cids_without_pKa]))
                continue
            
            data['dG0_r_tag'] = np.vstack([data['dG0_r_tag'], nist_row_data.dG0_r])
            data['pH'] = np.vstack([data['pH'], nist_row_data.pH])
            data['I'] = np.vstack([data['I'], nist_row_data.I])
            data['pMg'] = np.vstack([data['pMg'], nist_row_data.pMg])
            data['T'] = np.vstack([data['T'], nist_row_data.T])
            ddG = self.ReverseTransformReaction(nist_row_data.reaction, 
                nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
                nist_row_data.T)
            data['ddG0_r'] = np.vstack([data['ddG0_r'], ddG])
            
            stoichiometric_row = np.zeros((1, len(data['cids_to_estimate'])))
            for cid, coeff in nist_row_data.reaction.iteritems():
                stoichiometric_row[0, data['cids_to_estimate'].index(cid)] = coeff
            
            data['S'] = np.vstack([data['S'], stoichiometric_row])
        
        data['dG0_r'] = data['dG0_r_tag'] - data['ddG0_r']
        
        # remove the columns that are all-zeros in S
        nonzero_columns = np.sum(abs(data['S']), 0).nonzero()[0]
        data['S'] = data['S'][:, nonzero_columns]
        data['cids_to_estimate'] = np.array(data['cids_to_estimate'])
        data['cids_to_estimate'] = data['cids_to_estimate'][nonzero_columns]
        
        return data
    
    def ReverseTransformReaction(self, reaction, pH, I, pMg, T):
        """
            Calculates the difference between dG'0_r and dG0_r
        """
        return sum([coeff * self.ReverseTransformCompound(cid, pH, I, pMg, T) \
                    for cid, coeff in reaction.iteritems()])
        
    def ReverseTransformCompound(self, cid, pH, I, pMg, T):
        """
            Calculates the difference between dG'0_f and dG0_f
        """
        return self.Transform(cid, pH, I, pMg, T)
    
    def GetPseudoisomerMap(self, cid):
        return self.GetDissociationTable(cid).GetPseudoisomerMap()
    
    def ConvertPseudoisomer(self, cid, dG0, nH_from, 
                            nH_to=None, nMg_from=0, nMg_to=0):
        #try:
        diss = self.GetDissociationTable(cid)
        if nH_to == None:
            nH_to = diss.min_nH
        return diss.ConvertPseudoisomer(dG0, nH_from, nH_to, nMg_from, nMg_to)
        #except KeyError as e:
        #    raise KeyError("Cannot find the pKas of C%05d (%s): nH_from = %s, nH_to = %s" % \
        #                   (cid, self.kegg.cid2name(cid), str(nH_from), str(nH_to)))

        
    def Transform(self, cid, pH, I, pMg, T):
        return self.GetDissociationTable(cid).Transform(pH, I, pMg, T)

###############################################################################

class DissociationTable(object):
    
    def __init__(self, cid=None):
        # ddGs is a dictionary whose keys are 4-tuples of (nH_above, nH_below, nMg_above, nMg_below)
        # and the values are pairs of (ddG0, reference)
        self.ddGs = {}
        
        # smiles_dict has the same keys are ddGs, but the values are pairs of
        # (smiles_above, smiles_below)
        self.smiles_dict = {}
        
        self.cid = cid
        self.min_nH = None # the nH of the most basic pseudoisomer
        self.min_charge = None # the charge of the most basic pseudoisomer
        self.min_dG0 = 0 # the dG0 of the most basic pseudoisomer
        self.kegg = Kegg.getInstance()

    def __str__(self):
        s = "Base: nH=%d, z=%d, dG0=%.1f kJ/mol\n" % \
            (self.min_nH, self.min_charge, self.min_dG0)
        for (nH_above, nH_below, nMg_above, nMg_below), (ddG, ref) in self.ddGs.iteritems():
            if nH_above != nH_below:
                s += "nH (%2d -> %2d) : %.1f kJ/mol [%s]\n" % (nH_above, nH_below, ddG, ref)
            if nMg_above != nMg_below:
                s += "nMg (%2d -> %2d) : %.1f kJ/mol [%s]\n" % (nMg_above, nMg_below, ddG, ref)
        return s

    def __iter__(self):
        return self.ddGs.__iter__()
    
    def ToDatabase(self, db, table_name):
        """
            Write to a table with this structure:
                cid INT, name TEXT, 
                nH_below INT, nH_above INT, 
                nMg_below INT, nMg_above INT, 
                smiles_below TEXT, smiles_above TEXT, 
                ddG REAL, ref TEXT
        """
        name = self.kegg.cid2name(self.cid)
        for key, (ddG, ref) in self.ddGs.iteritems():
            (nH_above, nH_below, nMg_above, nMg_below) = key
            smiles_below, smiles_above = self.smiles_dict.get(key, (None, None))
            db.Insert(table_name, [self.cid, name, nH_above, nH_below, 
                nMg_above, nMg_below, smiles_below, smiles_above, ddG, ref])

    def AddpKa(self, pKa, nH_below, nH_above, nMg=0, ref="", T=default_T, 
               smiles_below=None, smiles_above=None):
        if nH_below != nH_above+1:
            raise Exception('A H+ dissociation constant (pKa) has to represent an '
                            'increase of exactly one hydrogen: nH_below = %d, nH_above = %d' %
                            (nH_below, nH_above))
        
        ddG0 = R * T * np.log(10) * pKa
        key = (nH_above, nH_below, nMg, nMg)
        self.ddGs[key] = (-ddG0, ref) # adding H+ decreases dG0
        self.smiles_dict[key] = (smiles_above, smiles_below)
        
    def AddpKMg(self, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T, 
                smiles_below=None, smiles_above=None):
        if nMg_below != nMg_above+1:
            raise Exception('A Mg+2 dissociation constant (pK_Mg) has to represent an '
                            'increase of exactly one magnesium ion: nMg_below = %d, nMg_above = %d' %
                            (nMg_below, nMg_above))
        ddG0 = R * T * np.log(10) * pKMg - dG0_f_Mg
        key = (nH, nH, nMg_above, nMg_below)
        self.ddGs[key] = (-ddG0, ref) # adding Mg+2 decreases dG0
        self.smiles_dict[key] = (smiles_above, smiles_below)
    
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
            raise KeyError('The dissociation constant for C%05d: (nH=%d,nMg=%d) -> '
                            '(nH=%d,nMg=%d) is missing' % (self.cid, nH_from, nMg_from, nH_to, nMg_to))

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
    
    def SetFormationEnergyByNumHydrogens(self, dG0, nH):
        """ Uses the value of any pseudoisomer to set the base value of dG0 """
        self.min_dG0 = self.ConvertPseudoisomer(dG0, nH, self.min_nH)
        
    def SetFormationEnergyByCharge(self, dG0, charge):
        """ Uses the value of any pseudoisomer to set the base value of dG0 """
        nH = self.min_nH + (charge - self.min_charge)
        self.min_dG0 = self.ConvertPseudoisomer(dG0, nH, self.min_nH)
    
    def SetTransformedFormationEnergy(self, dG0_tag, pH, I, pMg, T):
        """ Sets the min_dG0 according to a transformed formation energy. """
        self.min_dG0 += dG0_tag - self.Transform(pH, I, pMg, T)
    
    def CalculateCharge(self):
        """ Calculate the charge for the most basic species """
        # get the charge and nH of the default pseudoisomer in KEGG:
        nH_z_pair = self.kegg.cid2nH_and_charge(self.cid)
        if nH_z_pair:
            nH, z = nH_z_pair
            self.min_charge = z + (self.min_nH - nH)
        else:
            self.min_charge = 0
        
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
                      nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + \
                      nH  * (R*T*np.log(10)*pH + DH) - (z**2) * DH
            dG0_tag_vec.append(dG0_tag)
        
        dG0_tag_total = -R * T * log_sum_exp([g / (-R*T) for g in dG0_tag_vec])
        
        return dG0_tag_total
        
    def GetPseudoisomerMap(self):
        pmap = PseudoisomerMap()
        for pdata in self.GenerateAll():
            pmap.Add(nH=pdata.hydrogens, z=pdata.net_charge, 
                     nMg=pdata.magnesiums, dG0=pdata.dG0)
        return pmap

###############################################################################

if (__name__ == '__main__'):
    db = SqliteDatabase("../res/gibbs.sqlite")
    dissociation = DissociationConstants.FromFile()
    dissociation.ToDatabase(db, 'dissociation_constants')
