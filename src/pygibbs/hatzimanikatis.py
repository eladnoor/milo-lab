from pylab import arange
from thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from thermodynamic_constants import J_per_cal
import pseudoisomer
import csv, sys
from pygibbs.pseudoisomers_data import DissociationTable
from pygibbs.kegg import Kegg
import logging
from pygibbs.thermodynamic_constants import R
import pylab
from toolbox.database import SqliteDatabase

HATZI_CSV_FNAME = "../data/thermodynamics/hatzimanikatis_cid.csv"

class Hatzi (Thermodynamics):
    
    def __init__(self, use_pKa=True):
        Thermodynamics.__init__(self)
        
        # the conditions in which Hatzimanikatis makes his predictions
        self.pH = 7.0
        self.I = 0.0
        self.pMg = 10.0
        self.T = 298.15
        
        self.kegg = Kegg.getInstance()
        
        # for some reason, Hatzimanikatis doesn't indicate that H+ is zero,
        # so we add it here
        H_pmap = pseudoisomer.PseudoisomerMap()
        H_pmap.Add(0, 0, 0, 0)
        self.cid2pmap_dict = {80 : H_pmap}
        
        self.cid2dG0_tag_dict = {80 : 0}
        self.cid2charge_dict = {80 : 0}
        
        if use_pKa:
            self.cid2DissociationTable = DissociationTable.ReadDissociationCsv()
        else:
            self.cid2DissociationTable = {}
        
        for row in csv.DictReader(open(HATZI_CSV_FNAME, 'r')):
            cid = int(row['ENTRY'][1:])
            self.cid2source_string[cid] = 'Jankowski et al. 2008'
            if row['DELTAG'] == "Not calculated":
                continue
            self.cid2dG0_tag_dict[cid] = float(row['DELTAG']) * J_per_cal
            self.cid2charge_dict[cid] = int(row['CHARGE'])

    def charge2nH(self, cid, charge):
        base_charge = self.kegg.cid2charge(cid) or 0
        base_nH = self.kegg.cid2num_hydrogens(cid) or 0
        return base_nH + (charge - base_charge)

    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        elif cid in self.cid2dG0_tag_dict:
            charge = self.cid2charge_dict[cid]
            try:
                nH = self.charge2nH(cid, charge)
            except KeyError:
                logging.warning('The compound C%05d is missing from KEGG' % cid)
                nH = 0
            
            # Unlike Alberty, Hatzimanikatis writes down the transformed 
            # formation energies but leaves the potential of H+ to be not 0.
            # We therefore need to subtract that potential from 
            # the given formation energies (we already set H+ to 0).
            ddG0 = nH * R * self.T * pylab.log(10) * self.pH
            dG0_tag = self.cid2dG0_tag_dict[cid] + ddG0
            
            if cid not in self.cid2DissociationTable:
                diss_table = DissociationTable(cid)
                diss_table.min_nH = nH
                diss_table.min_charge = charge
                self.cid2DissociationTable[cid] = diss_table

            self.cid2DissociationTable[cid].SetTransformedFormationEnergy(
                dG0_tag, pH=7.0, I=0.0, pMg=10.0, T=298.15)
            self.cid2pmap_dict[cid] = self.cid2DissociationTable[cid].GetPseudoisomerMap()
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy(
                "The compound C%05d does not have a value for its"
                " formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2dG0_tag_dict.keys())
        
if (__name__ == "__main__"):
    db = SqliteDatabase('../res/gibbs.sqlite')
    H = Hatzi()
    #H.ToDatabase(db, 'hatzi_gc')
    #H.I = 0.25
    #H.T = 300;
    #sparse_reaction = {13:-1, 1:-1, 9:2}
    sparse_reaction = {36:-1, 3981:1}
    #sys.stdout.write("The dG0_r of PPi + H20 <=> 2 Pi: \n\n")
    sys.stdout.write("%5s | %5s | %6s | %6s\n" % ("pH", "I", "T", "dG0_r"))
    for pH in arange(5, 9.01, 0.25):
        H.pH = pH
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % 
                         (H.pH, H.I, H.T, H.reaction_to_dG0(sparse_reaction)))

    print '-'*50
    print H.cid2DissociationTable[36]
    print H.cid2DissociationTable[3981]
    print H.cid2pmap_dict[36].Transform(pH=7, I=0, pMg=10, T=298.15)
    print H.cid2pmap_dict[3981].Transform(pH=7, I=0, pMg=10, T=298.15)
