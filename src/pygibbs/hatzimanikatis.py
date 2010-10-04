import csv, sys
from common import *
from pylab import arange

class Hatzi:
    def __init__(self):
        csv_reader = csv.reader(open("../data/hatzimanikatis_cid.csv", "r"))
        csv_reader.next()
        self.cid2dG0_f = {}
        self.cid2charge = {}
        for row in csv_reader:
            cid = int(row[0][1:])
            if (row[1] == "Not calculated"):
                continue
            self.cid2dG0_f[cid] = float(row[1]) * 4.184
            self.cid2charge[cid] = int(row[3])
            
    def reaction_to_dG0(self, sparse_reaction, pH=default_pH, I=default_I, T=default_T):
        """
            sparse_reaction must be a dictionary where the keys are CIDs (integers) and the values are the stoichimetric coeff.
            Substrates have negative values and products have positive values.

            Note: Hatzimanikatis does not provide the nH of each compound, but since there is only one pseudoisomer
                  for each compound, nH can be normalized to the neutral form (i.e. nH = 0 when the charge is 0).
                  And in this case nH will be equal to the charge.
        """
        dG0_r = 0
        charge_diff = 0
        for (cid, coeff) in sparse_reaction.iteritems():
            if (cid not in self.cid2dG0_f):
                raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)
            dG0_f = self.cid2dG0_f[cid]
            z = self.cid2charge[cid]
            nH = z
            dG0_r += coeff * transform(dG0_f, nH, z, pH, I, T) 
        return dG0_r
    
    def display_pmap(self, cid):
        print "C%05d | %2d | %6.2f" % (cid, self.cid2charge[cid], self.cid2dG0_f[cid])
            
if (__name__ == "__main__"):
    H = Hatzi()
    I = 0.25; T = 300;
    sparse_reaction = {13:-1, 1:-1, 9:2}
    sys.stdout.write("The dG0_r of PPi + H20 <=> 2 Pi: \n\n")
    sys.stdout.write("%5s | %5s | %6s | %6s\n" % ("pH", "I", "T", "dG0_r"))
    for pH in arange(5, 9.01, 0.25):
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % (pH, I, T, H.reaction_to_dG0(sparse_reaction, pH, I=I, T=T)))
