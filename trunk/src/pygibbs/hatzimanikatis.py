import csv, sys
from pylab import arange
from thermodynamics import Thermodynamics, MissingCompoundFormationEnergy, J_per_cal

class Hatzi (Thermodynamics):
    def __init__(self):
        Thermodynamics.__init__(self)
        csv_reader = csv.reader(open("../data/thermodynamics/hatzimanikatis_cid.csv", "r"))
        csv_reader.next()
        self.cid2pmap_dict = {80 : {(0, 0) : 0} } # for some reason, Hatzimanikatis doesn't indicate that H+ is zero
        for row in csv_reader:
            cid = int(row[0][1:])
            if (row[1] == "Not calculated"):
                continue
            dG0_f = float(row[1]) * J_per_cal
            z = int(row[3])
            # Hatzimanikatis does not give the number of hydrogen atoms in his table, so we use the charge instead.
            # This should work since the transformed energies are linear in nH, and therefore it shouldn't matter
            # when calculating dG0 or reactions (since the real number of hydrogens will cancel out).
            nH = z
            self.cid2pmap_dict[cid] = {(nH, z) : dG0_f}
   
    def cid2pmap(self, cid):
        if (cid in self.cid2pmap_dict):
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

        
if (__name__ == "__main__"):
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
