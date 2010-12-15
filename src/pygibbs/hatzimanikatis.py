from pylab import arange
from thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from thermodynamic_constants import J_per_cal
import pseudoisomer
from toolbox.util import ReadCsvWithTitles

class Hatzi (Thermodynamics):
    def __init__(self):
        Thermodynamics.__init__(self)
        self.source_string = "Hatzimanikatis"
        H_pmap = pseudoisomer.PseudoisomerMap()
        H_pmap.Add(0, 0, 0, 0)
        self.cid2pmap_dict = {80 : H_pmap} # for some reason, Hatzimanikatis doesn't indicate that H+ is zero
        for row in ReadCsvWithTitles("../data/thermodynamics/hatzimanikatis_cid.csv"):
            cid = int(row['ENTRY'][1:])
            if (row['DELTAG'] == "Not calculated"):
                continue
            dG0_f = float(row['DELTAG']) * J_per_cal
            z = int(row['CHARGE'])
            # Hatzimanikatis does not give the number of hydrogen atoms in his table, so we use the charge instead.
            # This should work since the transformed energies are linear in nH, and therefore it shouldn't matter
            # when calculating dG0 or reactions (since the real number of hydrogens will cancel out).
            nH = z
            nMg = 0
            self.cid2pmap_dict[cid] = pseudoisomer.PseudoisomerMap()
            self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0_f)
   
    def cid2pmap(self, cid):
        if (cid in self.cid2pmap_dict):
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)

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
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % (H.pH, H.I, H.T, H.reaction_to_dG0(sparse_reaction)))

    sparse_reaction = {13:-1, 1:-1, 9:2}
