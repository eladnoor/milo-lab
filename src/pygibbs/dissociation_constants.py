import csv, logging, sys, time, threading
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamic_constants import R, default_T, dG0_f_Mg, debye_huckel,\
    default_pH, RedoxCarriers
import numpy as np
from toolbox.util import log_sum_exp
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.pseudoisomers_data import PseudoisomerEntry
from pygibbs.kegg_errors import KeggParseException
from optparse import OptionParser
from toolbox.molecule import OpenBabelError
from pygibbs.nist import Nist
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics

class MissingDissociationConstantError(Exception):
    
    def __init__(self, value, cid):
        self.value = value
        self.cid = cid
    def __str__(self):
        return repr(self.value)
     
class DissociationConstants(object):
    
    def __init__(self):
        self.cid2DissociationTable = {}
        self.kegg = Kegg.getInstance()
    
    @staticmethod
    def _CreateDatabase(db, table_name, drop_if_exists=True):
        db.CreateTable(table_name, """
            cid INT, name TEXT, 
            nH_below INT, nH_above INT, 
            nMg_below INT, nMg_above INT,
            charge_below INT, charge_above INT,
            mol_below TEXT, mol_above TEXT, 
            ddG REAL, ref TEXT""", drop_if_exists=drop_if_exists)
        
    
    def ToDatabase(self, db, table_name):
        """
            Load the data regarding pKa values according to KEGG compound IDs.
        """
        kegg = Kegg.getInstance()
        DissociationConstants._CreateDatabase(db, table_name)
        
        for cid in sorted(self.cid2DissociationTable.keys()):
            name = kegg.cid2name(cid)
            diss_table = self.cid2DissociationTable[cid]
            if diss_table is None:
                db.Insert(table_name, [cid, name] + [None] * 8)
            else:
                for row in diss_table.ToDatabaseRow():
                    db.Insert(table_name, [cid, name] + row)

        db.Commit()
            
    @staticmethod
    def FromFileToDB(file_name, db, table_name):
        """
            Parses a CSV file that contains pKa and pKMg data for many compounds
            and returns a dictionary of their DissociationTables, where the key
            is the CID.
            
            We support to CSV formats (legacy issues, sorry):
            1) cid, name, nH_below, nH_above, nMg_below, nMg_above, mol_below, mol_above, ddG, ref
            2) cid, name, type, T, nH_below, nH_above, nMg_below, nMg_above, mol_below, mol_above, pK, ref
        """
        
        kegg = Kegg.getInstance()
        DissociationConstants._CreateDatabase(db, table_name)

        for i, row in enumerate(csv.DictReader(open(file_name, 'r'))):
            if 'pK' not in row and 'ddG' not in row:
                raise Exception("The CSV file is not in a recognized format: "
                                "there should be a column named ddG or pK")
            try:
                if not row['cid']:
                    continue # without a CID we cannot match this to the dG0 table
                cid = int(row['cid'])
                name = row['name'] or kegg.cid2name(cid)
                logging.debug("Parsing row #%d, compound %s (C%05d)" %
                              (i, name, cid))
    
                nH_below = int(row['nH_below'])
                nH_above = int(row['nH_above'])
                nMg_below = int(row['nMg_below'])
                nMg_above = int(row['nMg_above'])
                mol_below = row['mol_below'] or None
                mol_above = row['mol_above'] or None
                ref = row['ref']
                
                if 'ddG' in row: # this is the 1st format type
                    ddG = float(row['ddG'])
                elif 'pK' in row: # this is the 2nd format type
                    pK = float(row['pK'] or 0)
                    T = float(row['T'] or default_T)
                    if row['type'] == 'acid-base':
                        if nMg_below != nMg_above or nH_below != nH_above+1:
                            raise Exception('wrong nMg and nH values')
                        ddG = -R * T * np.log(10) * pK
                    elif row['type'] == 'Mg':
                        if nMg_below != nMg_above+1 or nH_below != nH_above:
                            raise Exception('wrong nMg and nH values')
                        ddG = -R * T * np.log(10) * pK + dG0_f_Mg
                    elif row['type'] == '':
                        if nMg_below != nMg_above or nH_below != nH_above:
                            raise Exception('wrong nMg and nH values')
                        ddG = None
                    else:
                        raise Exception('unknown dissociation type: ' + row['type'])

            except Exception as e:
                raise Exception("Row %i: %s" % (i, str(e)))

            db.Insert(table_name, [cid, name, nH_below, nH_above, 
                                   nMg_below, nMg_above, mol_below,
                                   mol_above, ddG, ref])
        
        db.Commit()
    
    @staticmethod
    def FromChemAxon(cid2mol=None, html_writer=None):
        kegg = Kegg.getInstance()
        diss = DissociationConstants()
        if cid2mol is None:
            cid2mol = dict([(cid, None) for cid in kegg.get_all_cids()])
        
        for cid, mol in sorted(cid2mol.iteritems()):
            logging.info("Using ChemAxon to find the pKa values for %s - C%05d" %
                         (kegg.cid2name(cid), cid))
            if html_writer:
                html_writer.write('<h2>%s - C%05d</h2>\n' %
                                  (kegg.cid2name(cid), cid))
            # if this CID is not assigned to a Molecule, use the KEGG database
            # to create a Molecule for it.
            if mol is None:
                try:
                    mol = kegg.cid2mol(cid)
                except KeggParseException:
                    continue

            diss_table = mol.GetDissociationTable()
            diss.cid2DissociationTable[cid] = diss_table
            if diss_table and html_writer:
                diss_table.WriteToHTML(html_writer)
                html_writer.write('</br>\n')
        return diss
    
    @staticmethod
    def FromDatabase(db, table_name='dissociation_constants'):
        diss = DissociationConstants()

        for row in db.DictReader(table_name):
            cid = row['cid']
            if row['nH_below'] is None:
                diss.cid2DissociationTable[cid] = None
            else:
                if cid not in diss.cid2DissociationTable:
                    diss.cid2DissociationTable[cid] = DissociationTable(cid)
                try:
                    diss.cid2DissociationTable[cid].UpdateDatabaseRow(row)
                except OpenBabelError as e:
                    raise Exception("Cannot read this row from the database: "
                                    + str(row) + '\n' + str(e))
                    
        #diss.CalculateAllCharges()
        return diss
    
    @staticmethod
    def FromPublicDB():
        db = SqliteDatabase('../data/public_data.sqlite')
        return DissociationConstants.FromDatabase(db)
    
    def AddpKa(self, cid, pKa, nH_below, nH_above, nMg=0, ref="", T=default_T, 
               mol_below=None, mol_above=None):
        diss_table = self.GetDissociationTable(cid)
        diss_table.AddpKa(pKa, nH_below, nH_above, nMg,
                          ref, T, mol_below, mol_above)
        
    def AddpKMg(self, cid, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T, 
                mol_below=None, mol_above=None):
        diss_table = self.GetDissociationTable(cid)
        diss_table.AddpKMg(pKMg, nMg_below, nMg_above, nH, 
                          ref, T, mol_below, mol_above)
        
    def SetOnlyPseudoisomerMolecule(self, cid, mol, nMg=0):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            diss_table = DissociationTable(cid)
            self.cid2DissociationTable[cid] = diss_table
        diss_table.SetOnlyPseudoisomerMolecule(mol, nMg=nMg)

    def SetOnlyPseudoisomer(self, cid, nH, z, nMg=0):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            diss_table = DissociationTable(cid)
            self.cid2DissociationTable[cid] = diss_table
        diss_table.SetOnlyPseudoisomer(nH, z, nMg=nMg)
    
    def UpdateMinNumHydrogens(self, cid, min_nH):
        """
            Sets the minimal number of hydrogen atoms for a specified CID.
            
            If the value is not lower than the previously provided min_nH, the
            lower value is kept unchanged. 
        """
        diss_table = self.GetDissociationTable(cid)
        diss_table.UpdateMinNumHydrogens(min_nH)
            
    def GetDissociationTable(self, cid, create_if_missing=True):
        if cid not in self.cid2DissociationTable and create_if_missing:
            try:
                kegg = Kegg.getInstance()
                mol = kegg.cid2mol(cid)
                diss_table = DissociationTable.FromMolecule(mol)
            except KeggParseException:
                diss_table = None
            self.cid2DissociationTable[cid] = diss_table 
        
        return self.cid2DissociationTable.get(cid, None)
        
    def GetMol(self, cid, nH, nMg=0):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            return None
        else:
            return diss_table.GetMol(nH, nMg)
    
    def CalculateAllCharges(self):
        for diss_table in self.cid2DissociationTable.values():
            if diss_table is not None:
                try:
                    diss_table.CalculateCharge()
                except OpenBabelError as e:
                    raise Exception("error calculating charge for C%05d" % diss_table.cid
                                    + "\n" + str(e))

    def GetAllCids(self):
        return set(self.cid2DissociationTable.keys())

    def GetCid2nH_nMg(self, pH, I, pMg, T):
        cid2nH_nMg = {} # the nH that is to be used in the reverse transform
        for cid in self.GetAllCids():
            diss_table = self.GetDissociationTable(cid)
            if diss_table is not None:
                nH, nMg = diss_table.GetMostAbundantPseudoisomer(pH=pH, I=I, pMg=pMg, T=T)
                cid2nH_nMg[cid] = (nH, nMg)
            else:
                # assume nH=0 and nMg=0 by default for compounds without explicit formulas
                cid2nH_nMg[cid] = (0, 0)
        return cid2nH_nMg

    def ReverseTransformNistRows(self, nist_rows, cid2nH_nMg=None):
        all_cids = set()
        for nist_row_data in nist_rows:
            all_cids.update(nist_row_data.GetAllCids())
        all_cids = sorted(all_cids)
        
        data = {}
        data['dG0_r_tag'] = [] # the transformed free energy of the reactions dG'0_r
        data['dG0_r'] = [] # the chemical free energy of the reactions dG0_r
        data['ddG0_r'] = [] # dG'0_r - dG0_r  (which is only a function of the conditions and pKas)
        data['pH'] = []
        data['I'] = []
        data['pMg'] = []
        data['T'] = []
        data['S'] = np.zeros((len(all_cids), 0)) # stoichiometric matrix
        data['nist_rows'] = [] # The NIST rows that were used in S (since some might
                               # be dropped in the process, e.g. missing pKa)
        
        for nist_row_data in nist_rows:
            # check that all participating compounds have a known pKa
            dG0_prime = nist_row_data.dG0_r
            try:
                dG0 = self.ReverseTransformNistRow(nist_row_data, cid2nH_nMg)
            except MissingDissociationConstantError as e:
                logging.debug('Reaction %s involves compounds (C%05d) with missing pKa '
                              'values: %s' % (nist_row_data.ref_id, 
                                              e.cid, str(nist_row_data.reaction)))
                continue
            
            data['dG0_r_tag'].append(dG0_prime)
            data['pH'].append(nist_row_data.pH)
            data['I'].append(nist_row_data.I)
            data['pMg'].append(nist_row_data.pMg)
            data['T'].append(nist_row_data.T)
            data['nist_rows'].append(nist_row_data)
            data['ddG0_r'].append(dG0_prime - dG0)
            data['dG0_r'].append(dG0)
            
            # convert the reaction's sparse representation to a row vector
            stoichiometric_row = np.zeros((len(all_cids), 1))
            for cid, coeff in nist_row_data.reaction.iteritems():
                stoichiometric_row[all_cids.index(cid), 0] = coeff
            data['S'] = np.hstack([data['S'], stoichiometric_row])
        
        # remove the compounds that are all-zeros in S
        nonzero_rows = abs(data['S']).sum(1).nonzero()[0]
        data['S'] = np.matrix(data['S'][nonzero_rows, :])
        data['cids_to_estimate'] = [all_cids[i] for i in nonzero_rows]
        
        return data
    
    def ReverseTransformNistRow(self, nist_row_data, cid2nH_nMg=None, suppress_missing_pka_exception=False):
        """
            Given a NistRowData object, returns the reverse Legendre transform
            value for its dG0.
        """
        ddG = self.ReverseTransformReaction(nist_row_data.reaction, 
            nist_row_data.pH, nist_row_data.I, nist_row_data.pMg,
            nist_row_data.T, cid2nH_nMg=cid2nH_nMg,
            suppress_missing_pka_exception=suppress_missing_pka_exception)
        return nist_row_data.dG0_r - ddG

    def ReverseTransformReaction(self, reaction, pH, I, pMg, T, cid2nH_nMg=None, suppress_missing_pka_exception=False):
        """
            Calculates the difference between dG0_prime and dG0
        """
        ddG0 = 0
        for cid, coeff in reaction.iteritems():
            diss_table = self.GetDissociationTable(cid)
            if diss_table is None:
                if suppress_missing_pka_exception:
                    continue
                # probably a compound without an implicit formula
                raise MissingDissociationConstantError("", cid)
            elif not cid2nH_nMg:
                ddG0 += coeff * diss_table.GetDeltaDeltaG0(pH, I, pMg, T)
            else:
                nH, nMg = cid2nH_nMg[cid]
                ddG0 += coeff * diss_table.GetDeltaDeltaG0(pH, I, pMg, T, nH=nH, nMg=nMg)
        
        return ddG0
    
    def GetPseudoisomerMap(self, cid):
        return self.GetDissociationTable(cid).GetPseudoisomerMap()
    
    def ConvertPseudoisomer(self, cid, dG0, nH_from, 
                            nH_to=None, nMg_from=0, nMg_to=0):
        diss = self.GetDissociationTable(cid)
        if nH_to == None:
            nH_to = diss.min_nH
        return diss.ConvertPseudoisomer(dG0, nH_from, nH_to, nMg_from, nMg_to)
        
    def Transform(self, cid, pH, I, pMg, T):
        return self.GetDissociationTable(cid).Transform(pH, I, pMg, T)

    def GetMostAbundantMol(self, cid, pH, I, pMg, T):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            return None
        return diss_table.GetMostAbundantMol(pH, I, pMg, T)
    
    def GetMostAbundantPseudoisomer(self, cid, pH, I, pMg, T):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            return None
        return diss_table.GetMostAbundantPseudoisomer(pH, I, pMg, T)

    def GetAnyMol(self, cid):
        diss_table = self.GetDissociationTable(cid)
        if diss_table is None:
            return None
        return diss_table.GetAnyMol()
        
    def WriteToHTML(self, html_writer):
        for cid in sorted(self.GetAllCids()):
            diss_table = self.GetDissociationTable(cid)
            html_writer.write('<h2>%s - C%05d</h2>\n' % (self.kegg.cid2name(cid), cid))
            if diss_table is not None:
                diss_table.WriteToHTML(html_writer)
                html_writer.write('</br>\n')

###############################################################################


class DissociationTable(object):
    
    def __init__(self, cid=None):
        # ddGs is a dictionary whose keys are 4-tuples of (nH_above, nH_below, nMg_above, nMg_below)
        # and the values are pairs of (ddG0, reference)
        self.ddGs = {}
        
        # mol_dict is a dictionary from pairs of (nH, nMg) to a pair
        # (smiles, mol) where "smiles" is a string and "mol" is a Molecule object.
        # the usage of "mol" is lazy, i.e. it will be created only if requested
        # specifically, otherwise it will be None
        self.mol_dict = {}
        
        self.cid = cid
        self.min_nH = None # the nH of the most basic pseudoisomer
        self.min_charge = None # the charge of the most basic pseudoisomer
        self.min_dG0 = 0 # the dG0 of the most basic pseudoisomer

    def __len__(self):
        return len(self.ddGs)

    def __str__(self):
        T = default_T
        s = "Base   nH=%d nMg=0 z=%d : dG0=%.1f kJ/mol\n" % \
            (self.min_nH, self.min_charge, self.min_dG0)
        for key in sorted(self.ddGs.keys()):
            nH_above, nH_below, nMg_above, nMg_below = key
            ddG, _ref = self.ddGs[key]
            if nH_above != nH_below:
                pKa = -ddG / (R * T * np.log(10))
                s += "pKa (nH=%d -> nH=%d) : %.1f\n" % \
                    (nH_above, nH_below, pKa)
            elif nMg_above != nMg_below:
                pKMg = (-ddG + dG0_f_Mg) / (R * T * np.log(10))
                s += "pKMg (nMg=%d -> nMg=%d) : %.1f\n" % \
                    (nMg_above, nMg_below, pKMg)

        for nH, nMg in self.mol_dict.keys():
            s += "Pseudoisomer nH=%d nMg=%d : %s\n" % \
                 (nH, nMg, self.GetMolString(nH, nMg))
        
        return s

    @staticmethod
    def FromMolecule(mol, cid=None, nMg=0):
        diss_table = DissociationTable(cid)
        diss_table.SetOnlyPseudoisomerMolecule(mol, nMg=nMg)
        return diss_table

    def WriteToHTML(self, html_writer, T=default_T, draw_svg=True):
        dict_list = []
        if not self.ddGs:
            nH = self.min_nH
            nMg = 0
            ddG = 0.0
            ref = ""
            d = {'nH below':nH, 'nH above':nH,
                 'nMg below':nMg, 'nMg above':nMg,
                 'ddG0':'%.1f' % ddG, 'reference':ref}
            dict_list.append(d)
        else:
            for (nH_above, nH_below, nMg_above, nMg_below), (ddG, ref) in self.ddGs.iteritems():
                d = {'nH below':nH_below, 'nH above':nH_above,
                     'nMg below':nMg_below, 'nMg above':nMg_above,
                     'ddG0':'%.1f' % ddG, 'reference':ref}
                if nH_below == nH_above+1:
                    d['pK<sub>a</sub>'] = -ddG / (R * T * np.log(10))
                elif nMg_below == nMg_above+1:
                    d['pK<sub>Mg</sub>'] = (-ddG + dG0_f_Mg) / (R * T * np.log(10))
                dict_list.append(d)
                

        for d in dict_list:
            if draw_svg:
                d['species below'] = self.GetSVG(d['nH below'], d['nMg below'])
                d['species above'] = self.GetSVG(d['nH above'], d['nMg above'])
            else:
                d['species below'] = self.GetMolString(d['nH below'], d['nMg below'])
                d['species above'] = self.GetMolString(d['nH above'], d['nMg above'])
        
        dict_list.sort(key=lambda(k):(k['nH below'], k['nMg below']))
        html_writer.write_table(dict_list, headers=['nH below', 'nH above', 
            'nMg below', 'nMg above', 'species below', 'species above',
            'ddG0', 'pK<sub>a</sub>', 'pK<sub>Mg</sub>', 'reference'])        

    def __iter__(self):
        return self.ddGs.__iter__()
    
    def GetMol(self, nH=None, nMg=0):
        from toolbox.molecule import Molecule
        
        if nH is None:
            nH = self.min_nH
        if (nH, nMg) not in self.mol_dict:
            return None
        s, mol = self.mol_dict[nH, nMg]
        if mol is None:
            mol = Molecule.FromSmiles(s)
        self.mol_dict[nH, nMg] = (s, mol)
        return mol
    
    def GetSVG(self, nH=None, nMg=0):
        mol = self.GetMol(nH, nMg)
        if mol is not None:
            return mol.ToSVG()
        else:
            return ""
    
    def GetMolString(self, nH=None, nMg=0):
        if nH is None:
            nH = self.min_nH
        if (nH, nMg) not in self.mol_dict:
            return None
        s, _ = self.mol_dict[nH, nMg]
        return s
    
    def SetMolString(self, nH, nMg, s):
        self.UpdateMinNumHydrogens(nH)
        if s:
            self.mol_dict[nH, nMg] = (str(s), None)
            
    def SetMol(self, nH, nMg, mol):
        self.mol_dict[nH, nMg] = (mol.ToSmiles(), mol)
    
    def ToDatabaseRow(self):
        """
            Return:
                A list of rows to insert into the database
        """
        res = []
        if not self.ddGs:
            nH = self.min_nH
            z = self.min_charge
            nMg = 0
            mol_str = self.GetMolString(nH, nMg)
            ddG = 0.0
            ref = ""
            res.append([nH, nH, nMg, nMg, z, z, mol_str, mol_str, ddG, ref])
        else:
            for key in sorted(self.ddGs.keys()):
                nH_above, nH_below, nMg_above, nMg_below = key
                mol_below = self.GetMolString(nH_below, nMg_below)
                mol_above = self.GetMolString(nH_above, nMg_above)
                ddG, ref = self.ddGs[key]
                if self.min_charge is not None and self.min_nH is not None:
                    z_below = self.min_charge + (nH_below - self.min_nH) + 2*nMg_below
                    z_above = self.min_charge + (nH_above - self.min_nH) + 2*nMg_above
                res.append([nH_below, nH_above, nMg_below, nMg_above, 
                            z_below, z_above, mol_below, mol_above, ddG, ref])
        return res

    def UpdateDatabaseRow(self, row):
        (nH_below, nMg_below) = (row['nH_below'], row['nMg_below'])
        (nH_above, nMg_above) = (row['nH_above'], row['nMg_above'])
        self.SetMolString(nH_above, nMg_above, row['mol_above'])
        self.SetMolString(nH_below, nMg_below, row['mol_below'])
        if self.min_charge is None and row['charge_below'] is not None:
            self.min_charge = row['charge_below'] - (nH_below - self.min_nH) - 2*nMg_below

        if (nH_below, nMg_below) == (nH_above, nMg_above):
            return

        key = (nH_above, nH_below, nMg_above, nMg_below)
        self.ddGs[key] = (row['ddG'], row['ref'])

    def AddpKa(self, pKa, nH_below, nH_above, nMg=0, ref="", T=default_T, 
               mol_below=None, mol_above=None):
        if nH_below != nH_above+1:
            raise Exception('A H+ dissociation constant (pKa) has to represent an '
                            'increase of exactly one hydrogen: nH_below = %d, nH_above = %d' %
                            (nH_below, nH_above))
        
        ddG0 = R * T * np.log(10) * pKa
        key = (nH_above, nH_below, nMg, nMg)
        self.ddGs[key] = (-ddG0, ref) # adding H+ decreases dG0
        if mol_above is not None:
            self.SetMol(nH_above, nMg, mol_above)
        if mol_below is not None:
            self.SetMol(nH_below, nMg, mol_below)
        self.UpdateMinNumHydrogens(nH_above)
        
    def AddpKMg(self, pKMg, nMg_below, nMg_above, nH, ref="", T=default_T, 
                mol_below=None, mol_above=None):
        if nMg_below != nMg_above+1:
            raise Exception('A Mg+2 dissociation constant (pK_Mg) has to represent an '
                            'increase of exactly one magnesium ion: nMg_below = %d, nMg_above = %d' %
                            (nMg_below, nMg_above))
        ddG0 = R * T * np.log(10) * pKMg - dG0_f_Mg
        key = (nH, nH, nMg_above, nMg_below)
        self.ddGs[key] = (-ddG0, ref) # adding Mg+2 decreases dG0
        if mol_above is not None:
            self.SetMol(nH, nMg_above, mol_above)
        if mol_below is not None:
            self.SetMol(nH, nMg_below, mol_below)
        self.UpdateMinNumHydrogens(nH)
    
    def SetOnlyPseudoisomerMolecule(self, mol, nMg=0):
        """
            For compound which have no known pKa or pKMg, this method can be used
            to set the parameters of the only pseudoisomer.
        """
        if len(self.ddGs):
            raise ValueError("You tried to set the only-pseudoisomer of a compound that has pKas/pKMgs")
        nH, z = mol.GetHydrogensAndCharge()
        self.SetMol(nH, nMg, mol)
        self.min_nH = nH
        if z is not None:
            self.SetCharge(nH, z, nMg)
    
    def SetOnlyPseudoisomer(self, nH, z, nMg=0):
        """
            For compound which have no known pKa or pKMg, this method can be used
            to set the parameters of the only pseudoisomer.
        """
        if len(self.ddGs):
            raise ValueError("You tried to set the only-pseudoisomer of a compound that has pKas/pKMgs")
        self.min_nH = nH
        self.SetCharge(nH, z, nMg)
    
    def UpdateMinNumHydrogens(self, min_nH):
        if self.min_nH is None:
            self.min_nH = min_nH
        elif self.min_nH > min_nH: 
            self.min_charge -= (self.min_nH - min_nH)
            self.min_nH = min_nH

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
                '(nH=%d,nMg=%d) is missing' % (self.cid or 0, nH_from, nMg_from, nH_to, nMg_to),
                self.cid)

        raise Exception('A dissociation constant can either represent a'
                ' change in only one hydrogen or magnesium')
    
    def ConvertPseudoisomer(self, dG0, nH_from, nH_to=None, nMg_from=0, nMg_to=0):
        if nH_to is None:
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
        if comp.ref is not None:
            comp.ref += ';' + total_ref
        else:
            comp.ref = total_ref
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
        mol = self.GetAnyMol()
        if mol is None:
            # get the charge and nH of the default pseudoisomer in KEGG:
            kegg = Kegg.getInstance()
            nH_z_pair = kegg.cid2nH_and_charge(self.cid)
        else:
            nH_z_pair = mol.GetHydrogensAndCharge()
        
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
            elif nMg_below == nMg_above + 1: # this is a pK_Mg
                nH = nH_below # it doesn't matter which, since nH_below == nH_above
                nMg = nMg_below # since we always start from nMg=0 and go up
            else: # this is the only pseudoisomer and therefore nH and nMg all equal
                nH = nH_below
                nMg = nMg_below
                    
            pseudoisomers[nH, nMg] = self.ConvertPseudoisomerEntry(pdata, nH, nMg)
            pseudoisomers[nH, nMg].ref = ref
            
        return pseudoisomers.values()

    def GetTransformedDeltaGs(self, pH, I, pMg, T, nH=None, nMg=0):
        """
            Return:
                a list of the pseudoisomers and their transformed dG.
                each member of the list is a tuple: (nH, z, nMg, dG0')
            
            Note:
                Set the dG0 of one the pseudoisomer [nH, nMg] to 0.
        """
        if nH is None:
            nH = self.min_nH
        pdata = PseudoisomerEntry(net_charge=(self.min_charge + nH - self.min_nH),
                                  hydrogens=nH, magnesiums=nMg, smiles="", dG0=0)
        
        pseudoisomer_matrix = []
        for pseudoisomer in self.GenerateAllPseudoisomerEntries(pdata):
            ps_nH = pseudoisomer.hydrogens
            ps_nMg = pseudoisomer.magnesiums
            ps_z = pseudoisomer.net_charge
            ps_dG0 = pseudoisomer.dG0
            
            DH = debye_huckel(I)
            dG0_tag = ps_dG0 + \
                      ps_nMg * (R*T*np.log(10)*pMg - dG0_f_Mg) + \
                      ps_nH  * (R*T*np.log(10)*pH + DH) - (ps_z**2) * DH
            pseudoisomer_matrix.append((ps_nH, ps_z, ps_nMg, dG0_tag))
        
        pseudoisomer_matrix.sort(key=lambda(x):x[3])
        return pseudoisomer_matrix
     
    def GetDeltaDeltaG0(self, pH, I, pMg, T, nH=None, nMg=0):
        """
            Return:
                the transformed ddG0 = dG0'_f - dG0_f

            Note:
                assume that the dG0_f of one of the psuedoisomers
                (according to the given nH) is 0
        """
        if nH is None:
            nH = self.min_nH
        pseudoisomer_matrix = self.GetTransformedDeltaGs(pH, I, pMg, T, nH=nH, nMg=nMg)
        ddG0_f = -R * T * log_sum_exp([dG0_tag / (-R*T) for (_nH, _z, _nMg, dG0_tag) in pseudoisomer_matrix])
        return ddG0_f
    
    def GetMostAbundantPseudoisomer(self, pH, I, pMg, T):
        pseudoisomer_matrix = self.GetTransformedDeltaGs(pH, I, pMg, T)
        pseudoisomer_matrix.sort(key=lambda(x):x[3])
        nH, _z, nMg, _dG0_tag = pseudoisomer_matrix[0] # return the psuedoisomer with the smallest dG0_tag
        return (nH, nMg)
    
    def GetMostAbundantMol(self, pH, I, pMg, T):
        nH, nMg = self.GetMostAbundantPseudoisomer(pH, I, pMg, T)
        return self.GetMol(nH, nMg)
        
    def GetAnyMol(self):
        for nH, nMg in self.mol_dict.keys():
            mol = self.GetMol(nH, nMg)
            if mol is not None:
                return mol
        return None
    
    def Transform(self, pH, I, pMg, T):    
        return self.min_dG0 + self.GetDeltaDeltaG0(pH, I, pMg, T, nH=self.min_nH)
    
    def GetPseudoisomerMap(self, nH=None, nMg=0):
        """
            Generate all pseudoisomers around one specified pseudoisomer.
        """
        if nH is None:
            nH = self.min_nH
        pdata = PseudoisomerEntry(net_charge=self.min_charge, hydrogens=self.min_nH,
            magnesiums=0, smiles="", dG0=self.min_dG0, cid=self.cid)
        pdata = self.ConvertPseudoisomerEntry(pdata, nH_to=nH, nMg_to=nMg)
        pdata.ref = ''
        pmap = PseudoisomerMap()
        for pdata in self.GenerateAllPseudoisomerEntries(pdata):
            pmap.Add(nH=pdata.hydrogens, z=pdata.net_charge, 
                     nMg=pdata.magnesiums, dG0=pdata.dG0,
                     ref=pdata.ref)
        return pmap

###############################################################################

class DissociationThreads(threading.Thread):
    
    def __init__(self, group=None, target=None, name=None, args=(), kwargs={}):
        """
            args should contain (cid, smiles, semaphore)
        """
        threading.Thread.__init__(self, group, target, name, args, kwargs)
        self.cid, self.smiles, self.semaphore, self.db_lock, self.options = args
    
    def run(self):
        from toolbox.molecule import Molecule
        
        self.semaphore.acquire()
        
        start_time = time.time()

        logging.debug("SMILES: " + self.smiles)
        diss_table = Molecule._GetDissociationTable(self.smiles, fmt='smiles',
            mid_pH=default_pH, min_pKa=0, max_pKa=14, T=default_T)
        logging.debug("Min charge: %d" % diss_table.min_charge)
        logging.debug("Min nH: %d" % diss_table.min_nH)
        
        elapsed_time = time.time() - start_time
        self.db_lock.acquire()
        db = SqliteDatabase(self.options.db_file)
        kegg = Kegg.getInstance()
        name = kegg.cid2name(self.cid)
        
        if diss_table is not None:
            for row in diss_table.ToDatabaseRow():
                db.Insert(self.options.table_name, [self.cid, name] + row)
        else:
            db.Insert(self.options.table_name, [self.cid, name] + [None] * 10)
        del db
        self.db_lock.release()

        logging.info("Completed C%05d, elapsed time = %.1f sec" %
                     (self.cid, elapsed_time))

        self.semaphore.release()

###############################################################################

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-n", "--nist", action="store_true",
                          dest="nist",
                          default=False,
                          help="Calculate pKas only for the compounds in NIST"
                               " (otherwise use entire KEGG database)")
    opt_parser.add_option("-o", "--override", action="store_true",
                          dest="override_table",
                          default=False,
                          help="Drop the DB table and start from scratch")
    opt_parser.add_option("-p", "--threads", action="store", type="int",
                          dest="n_threads",
                          default=1,
                          help="Number of threads to use in parallel for calculating pKas")
    opt_parser.add_option("-d", "--database", action="store",
                          dest="db_file",
                          default="../data/public_data.sqlite",
                          help="The SQLite database to write to")
    opt_parser.add_option("-t", "--table_name", action="store",
                          dest="table_name",
                          default="dissociation_constants",
                          help="The name of the DB table for the results")
    return opt_parser


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    db = SqliteDatabase(options.db_file)
    kegg = Kegg.getInstance()
    
    if options.override_table:
        db.Execute("DROP TABLE IF EXISTS " + options.table_name)
    
    DissociationConstants._CreateDatabase(db, options.table_name, drop_if_exists=options.override_table)

    cids_to_calculate = set()
    if options.nist:
        cids_to_calculate.update(Nist().GetAllCids())
        cids_to_calculate.update(RedoxCarriers().GetAllCids())
        
        ptable = PsuedoisomerTableThermodynamics.FromCsvFile("../data/thermodynamics/formation_energies.csv")
        cids_to_calculate.update(ptable.get_all_cids())
    else:
        cids_to_calculate.update(kegg.get_all_cids())

    for row in db.Execute("SELECT distinct(cid) FROM %s" % options.table_name):
        if row[0] in cids_to_calculate:
            cids_to_calculate.remove(row[0])
    
    cid2smiles_and_mw = {}
    for cid in cids_to_calculate:
        # the compound CO is a special case where the conversion from InChI
        # to SMILES fails, so we add a specific override for it only
        if cid == 237:
            cid2smiles_and_mw[cid] = ("[C-]#[O+]", 28)
            continue
        
        try:
            comp = kegg.cid2compound(cid)
            mol = comp.GetMolecule()
            cid2smiles_and_mw[cid] = (mol.ToSmiles(), mol.GetExactMass())
        except KeggParseException:
            logging.debug("%s (C%05d) has no SMILES, skipping..." %
                          (kegg.cid2name(cid), cid))
        except OpenBabelError:
            logging.debug("%s (C%05d) cannot be converted to SMILES, skipping..." %
                          (kegg.cid2name(cid), cid))
        
    # Do not recalculate pKas for CIDs that are already in the database
    cids_to_calculate = cid2smiles_and_mw.keys()
    cids_to_calculate.sort(key=lambda(cid):(cid2smiles_and_mw[cid][1], cid))
    
    db_lock = threading.Lock()
    semaphore = threading.Semaphore(options.n_threads)
    for cid in cids_to_calculate:
        smiles, _ = cid2smiles_and_mw[cid]
        if not smiles:
            logging.info("The following compound is blacklisted: C%05d" % cid)
            continue

        thread = DissociationThreads(group=None, target=None, name=None,
                                     args=(cid, smiles, semaphore, db_lock, options), kwargs={})
        thread.start()
        
if __name__ == '__main__':
    logging.getLogger('').setLevel(logging.DEBUG)
    main()
