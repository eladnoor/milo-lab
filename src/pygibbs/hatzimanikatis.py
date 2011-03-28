from pylab import arange
from thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from thermodynamic_constants import J_per_cal
import pseudoisomer
import csv
from pygibbs.pseudoisomers_data import DissociationTable
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.kegg import Kegg
import logging

HATZI_CSV_FNAME = "../data/thermodynamics/hatzimanikatis_cid.csv"

class Hatzi (Thermodynamics):
    
    def __init__(self, use_pKa=True):
        Thermodynamics.__init__(self)
        
        # for some reason, Hatzimanikatis doesn't indicate that H+ is zero,
        # so we add it here
        H_pmap = pseudoisomer.PseudoisomerMap()
        H_pmap.Add(0, 0, 0, 0)
        self.cid2pmap_dict = {80 : H_pmap}
        
        if use_pKa:
            cid2DissociationTable = DissociationTable.ReadDissociationCsv()
        
        for row in csv.DictReader(open(HATZI_CSV_FNAME, 'r')):
            cid = int(row['ENTRY'][1:])

            self.cid2source_string[cid] = 'Jankowski et al. 2008'
            if row['DELTAG'] == "Not calculated":
                continue
            dG0_f = float(row['DELTAG']) * J_per_cal
            charge = int(row['CHARGE'])
            
            kegg = Kegg.getInstance()
            if cid not in cid2DissociationTable:
                self.cid2pmap_dict[cid] = PseudoisomerMap()
                try:
                    base_charge = kegg.cid2charge(cid) or 0
                    base_nH = kegg.cid2num_hydrogens(cid) or 0
                    nH = base_nH + (charge - base_charge)
                    self.cid2pmap_dict[cid].Add(nH=nH, z=charge, nMg=0, dG0=dG0_f)
                except KeyError:
                    logging.warning('the compound C%05d is missing from KEGG' % cid)
            else:
                diss_table = cid2DissociationTable[cid]
                diss_table.SetFormationEnergyByCharge(dG0_f, charge)
                self.cid2pmap_dict[cid] = diss_table.GetPseudoisomerMap()
   
    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy(
                "The compound C%05d does not have a value for its"
                " formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

        
if (__name__ == "__main__"):
    import sys

    H = Hatzi()
    H.I = 0.25
    H.T = 300;
    sparse_reaction = {13:-1, 1:-1, 9:2}
    sys.stdout.write("The dG0_r of PPi + H20 <=> 2 Pi: \n\n")
    sys.stdout.write("%5s | %5s | %6s | %6s\n" % ("pH", "I", "T", "dG0_r"))
    for pH in arange(5, 9.01, 0.25):
        H.pH = pH
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % 
                         (H.pH, H.I, H.T, H.reaction_to_dG0(sparse_reaction)))

    sparse_reaction = {13:-1, 1:-1, 9:2}
