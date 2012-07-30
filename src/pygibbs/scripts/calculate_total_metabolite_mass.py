import matplotlib.pyplot as plt
import numpy as np
import csv
from pygibbs.kegg import Kegg

def main():
    kegg = Kegg.getInstance()
    mass_conc = []
    names = []
    for row in csv.DictReader(open('../data/thermodynamics/abundance_bennett.csv', 'r')):
        if row['KEGG ID']:
            cid = int(row['KEGG ID'])
            conc = float(row['glucose'])
            MW = kegg.cid2compound(cid).mass
            if MW is not None:
                mass_conc += [conc * MW]
                names += [kegg.cid2name(cid)]
    #plt.bar(xrange(len(mass_conc)), mass_conc, 0.5)
    plt.pie(mass_conc, labels=names)
    plt.show()
    print np.sum(mass_conc)
        
        
        
if __name__ == "__main__":
    main()