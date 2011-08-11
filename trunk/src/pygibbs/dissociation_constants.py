import csv, logging
from kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamic_constants import R, default_T, dG0_f_Mg, debye_huckel
import numpy as np
from toolbox.util import log_sum_exp
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.pseudoisomers_data import PseudoisomerEntry

class MissingDissociationConstantError(Exception):
    pass

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
            try:
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
                else:
                    diss.SetOnlyPseudoisomer(cid, smiles_below, nH_below, nMg_below)
            except ValueError as e:
                raise ValueError("At row %i: %s" % (i, str(e)))
            except TypeError as e:
                raise TypeError("At row %i: %s" % (i, str(e)))
        
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
        
    def SetOnlyPseudoisomer(self, cid, smiles, nH, nMg=0):
        diss_table = self.GetDissociationTable(cid)
        diss_table.SetOnlyPseudoisomer(smiles, nH, nMg)
    
    def SetMinNumHydrogens(self, cid, min_nH):
        """
            Sets the minimal number of hydrogen atoms for a specified CID.
            
            If the value is not lower than the previously provided min_nH, the
            lower value is kept unchanged. 
        """
        diss_table = self.GetDissociationTable(cid)
        if not diss_table.min_nH or diss_table.min_nH > min_nH: 
            diss_table.min_nH = min_nH
            
    def GetDissociationTable(self, cid, create_if_missing=True):
        if cid not in self.cid2DissociationTable:
            if create_if_missing:
                diss = DissociationTable(cid)
                self.cid2DissociationTable[cid] = diss
                return diss
            else:
                return None
        else:
            return self.cid2DissociationTable[cid]
    
    def GetSmiles(self, cid, nH, nMg=0):
        if cid not in self.cid2DissociationTable:
            return None
        else:
            return self.GetDissociationTable(cid).GetSmiles(nH, nMg)
    
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

    def ReverseTranformNistRows(self, nist_rows, cid2nH=None, assume_no_pka_by_default=False):
        kegg = Kegg.getInstance()

        encountered_cids = self.GetAllCids()

        if assume_no_pka_by_default:
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
        data['nist_rows'] = np.zeros((0, 1)) # the index of the corresponding row in nist_rows
        
        for nist_row, nist_row_data in enumerate(nist_rows):
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
            data['nist_rows'] = np.vstack([data['nist_rows'], nist_row])
            ddG = self.ReverseTransformReaction(nist_row_data.reaction, 
                nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
                nist_row_data.T, cid2nH=cid2nH)
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
    
    def ReverseTransformReaction(self, reaction, pH, I, pMg, T, cid2nH=None):
        """
            Calculates the difference between dG'0_r and dG0_r
        """
        ddG0 = 0
        for cid, coeff in reaction.iteritems():
            diss = self.GetDissociationTable(cid)
            if not cid2nH:
                ddG0 += coeff * diss.GetDeltaDeltaG0(pH, I, pMg, T)
            else:
                ddG0 += coeff * diss.GetDeltaDeltaG0(pH, I, pMg, T, nH=cid2nH[cid])
        
        return ddG0
        
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
    
    def GetSmiles(self, nH=None, nMg=0):
        if nH == None:
            nH = self.min_nH
        return self.smiles_dict.get((nH, nMg), None)
    
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
            smiles_below = self.GetSmiles(nH_below, nMg_below)
            smiles_above = self.GetSmiles(nH_above, nMg_above)
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
        self.smiles_dict[nH_above, nMg] = smiles_above
        self.smiles_dict[nH_below, nMg] = smiles_below
        
        if self.min_nH is None or self.min_nH > nH_above:
            self.min_nH = nH_above
        
    def AddpKMg(self, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T, 
                smiles_below=None, smiles_above=None):
        if nMg_below != nMg_above+1:
            raise Exception('A Mg+2 dissociation constant (pK_Mg) has to represent an '
                            'increase of exactly one magnesium ion: nMg_below = %d, nMg_above = %d' %
                            (nMg_below, nMg_above))
        ddG0 = R * T * np.log(10) * pKMg - dG0_f_Mg
        key = (nH, nH, nMg_above, nMg_below)
        self.ddGs[key] = (-ddG0, ref) # adding Mg+2 decreases dG0
        self.smiles_dict[nH, nMg_above] = smiles_above
        self.smiles_dict[nH, nMg_below] = smiles_below

        if self.min_nH is None or self.min_nH > nH:
            self.min_nH = nH
    
    def SetOnlyPseudoisomer(self, smiles, nH, nMg=0):
        """
            For compound which have no known pKa or pKMg, this method can be used
            to set the parameters of the only pseudoisomer.
        """
        if len(self.ddGs):
            raise ValueError("You tried to set the only-pseudoisomer of a compound that has pKas/pKMgs")
        self.smiles_dict[nH, nMg] = smiles
        self.min_nH = nH
    
    def GetSingleStep(self, nH_from, nH_to, nMg_from, nMg_to):
        if nH_from == nH_to and nMg_from == nMg_to:
            return 0, None
        
        if nH_from != nH_to and nMg_from != nMg_to:
            raise Exception('A dissociation constant can either represent a'
                ' change in hydrogens or in magnesiums, but not both')
        
        try:
            if nMg_to == nMg_from+1 or nH_to == nH_from+1:
                ddG0, ref = self.ddGs[nH_from, nH_to, nMg_from, nMg_to]
                return ddG0, ref
            if nMg_from == nMg_to+1:
                ddG0, ref = self.ddGs[nH_from, nH_to, nMg_to, nMg_from]
                return -ddG0, ref
            if nH_from == nH_to+1:
                ddG0, ref = self.ddGs[nH_to, nH_from, nMg_from, nMg_to]
                return -ddG0, ref
        except KeyError:
            raise MissingDissociationConstantError(
                'The dissociation constant for C%05d: (nH=%d,nMg=%d) -> '
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
    
    def SetFormationEnergyByNumHydrogens(self, dG0, nH, nMg=0):
        """ Uses the value of any pseudoisomer to set the base value of dG0 """
        self.min_dG0 = self.ConvertPseudoisomer(dG0, nH_from=nH, 
            nH_to=self.min_nH, nMg_from=nMg, nMg_to=0)
        
    def SetFormationEnergyByCharge(self, dG0, charge, nMg=0):
        """ Uses the value of any pseudoisomer to set the base value of dG0 """
        nH = self.min_nH + (charge - self.min_charge)
        self.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
    
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
            
    def SetCharge(self, nH, z, nMg=0):
        self.min_charge = z + (self.min_nH - nH) - 2 * nMg
        
    def GenerateAll(self):
        if self.min_charge == None:
            raise Exception('The minimal charge has to be set before generating'
            ' all psuedoisomers')
        if self.min_dG0 == None:
            raise Exception('The base formation energy has to be set before generating'
            ' all psuedoisomers')
            
        pdata = PseudoisomerEntry(net_charge=self.min_charge, hydrogens=self.min_nH,
            magnesiums=0, smiles="", dG0=self.min_dG0, cid=self.cid)
        return self.GenerateAllPseudoisomerEntries(pdata)
    
    def GenerateAllPseudoisomerEntries(self, pdata):
        # Having nMg=0 is not necessary but makes everything easier,
        # therefore I assume it is.
        if pdata.magnesiums != 0:
            raise ValueError("Must start with nMg=0 when generating all pseudoisomers")
        
        pseudoisomers = {}
        pseudoisomers[pdata.hydrogens, pdata.magnesiums] = pdata.Clone()

        for (nH_above, nH_below, nMg_above, nMg_below), (_ddG0, ref) in self.ddGs.iteritems():
            if nH_below == nH_above + 1: # this is a pKa
                nMg = nMg_below # it doesn't matter which, since nMg_below == nMg_above
                if nH_below > pdata.hydrogens: 
                    nH = nH_below # creating a pseudoisomer with more nH than pdata
                elif nH_above < pdata.hydrogens:
                    nH = nH_above # creating a pseudoisomer with less nH than pdata
            elif nMg_below == nMg_above + 1:
                nH = nH_below # it doesn't matter which, since nH_below == nH_above
                nMg = nMg_below # since we always start from nMg=0 and go up
                    
            pseudoisomers[nH, nMg] = self.ConvertPseudoisomerEntry(pdata, nH, nMg)
            pseudoisomers[nH, nMg].ref = ref
            
        return pseudoisomers.values()

    def GetTransformedDeltaGs(self, pH, I, pMg, T, nH=None):
        """
            Return:
                a list of the pseudoisomers and their transformed dG.
                each member of the list is a tuple: (nH, z, nMg, dG0')
            
            Note:
                assume that the dG0_f of one of the psuedoisomers
                (according to the given nH) is 0
        """ 
        if nH is None:
            nH = self.min_nH
        z = self.min_charge + (nH - self.min_nH)
        
        pdata = PseudoisomerEntry(net_charge=z, hydrogens=nH, magnesiums=0, 
                                  smiles="", dG0=0)
        
        pseudoisomer_matrix = []
        for pseudoisomer in self.GenerateAllPseudoisomerEntries(pdata):
            nH = pseudoisomer.hydrogens
            nMg = pseudoisomer.magnesiums
            z = pseudoisomer.net_charge
            dG0 = pseudoisomer.dG0
            
            DH = debye_huckel(I)
            dG0_tag = dG0 + \
                      nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + \
                      nH  * (R*T*np.log(10)*pH + DH) - (z**2) * DH
            pseudoisomer_matrix.append((nH, z, nMg, dG0_tag))
        return pseudoisomer_matrix
     
    def GetDeltaDeltaG0(self, pH, I, pMg, T, nH=None):
        """
            Return:
                the transformed ddG0 = dG0'_f - dG0_f

            Note:
                assume that the dG0_f of one of the psuedoisomers
                (according to the given nH) is 0
        """ 
        pseudoisomer_matrix = self.GetTransformedDeltaGs(pH, I, pMg, T, nH=nH)
        ddG0_f = -R * T * log_sum_exp([dG0_tag / (-R*T) for (_nH, _z, _nMg, dG0_tag) in pseudoisomer_matrix])
        return ddG0_f
    
    def GetMostAbundantPseudoisomer(self, pH, I, pMg, T):
        pseudoisomer_matrix = self.GetTransformedDeltaGs(pH, I, pMg, T)
        pseudoisomer_matrix.sort(key=lambda(x):x[3])
        nH, _z, nMg, _dG0_tag = pseudoisomer_matrix[0] # return the psuedoisomer with the smallest dG0_tag
        return (nH, nMg)
    
    def Transform(self, pH, I, pMg, T):    
        return self.min_dG0 + self.GetDeltaDeltaG0(pH, I, pMg, T, nH=self.min_nH)
    
    def GetPseudoisomerMap(self):
        pmap = PseudoisomerMap()
        for pdata in self.GenerateAll():
            pmap.Add(nH=pdata.hydrogens, z=pdata.net_charge, 
                     nMg=pdata.magnesiums, dG0=pdata.dG0,
                     ref=pdata.ref)
        return pmap

###############################################################################

if __name__ == '__main__':
    db = SqliteDatabase("../res/gibbs.sqlite")
    dissociation = DissociationConstants.FromFile()
    dissociation.ToDatabase(db, 'dissociation_constants')
