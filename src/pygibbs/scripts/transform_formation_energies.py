from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
import csv
from pygibbs.kegg import Kegg

FormationEnergyFileName = "../data/thermodynamics/formation_energies.csv"

def main():
    ptable = PsuedoisomerTableThermodynamics.FromCsvFile(FormationEnergyFileName, label='testing')
    kegg = Kegg.getInstance()
    pH, I, pMg, T = (7.0, 0.25, 14, 298.15)
    
    output_csv = csv.writer(open('../res/formation_energies_transformed.csv', 'w'))
    output_csv.writerow(["cid","name","dG'0","pH","I","pMg","T",
                         "anchor","compound_ref","remark"])
    for cid in ptable.get_all_cids():
        pmap = ptable.cid2PseudoisomerMap(cid)
        dG0_prime = pmap.Transform(pH=pH, I=I, pMg=pMg, T=T)
        output_csv.writerow([cid, kegg.cid2name(cid), "%.1f" % dG0_prime, pH, I, pMg, T,
                             1, ptable.cid2source_string[cid]])


if __name__ == "__main__":
    main()