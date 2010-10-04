import csv, re, sys
from pylab import arange
from common import *

class Alberty:
    def read_alberty_mathematics(self, fname):
        """
            example line:
            acetatesp={{-369.31,-486.01,-1,3},{-396.45,-485.76,0,4}};
            
            the order of values is: (dG0, dH0, z, nH)
        """
        file = open(fname, 'r')
        alberty_name_to_pmap = {}
        for line in file.readlines():
            line.rstrip()
            if (line.find('=') == -1):
                continue
            (alberty_name, values) = line.split('sp=', 1)
            pmap = {}
            for token in re.findall("{([0-9\-\.\,_\s]+)}", values):
                val_list = token.split(',', 3)
                dG0 = float(val_list[0])
                z = int(val_list[2])
                nH = int(val_list[3])
                if (alberty_name.find("coA") != -1):
                    nH += 35
                pmap[(nH, z)] = dG0
            alberty_name_to_pmap[alberty_name] = pmap
        return alberty_name_to_pmap
        
        
    def read_alberty_kegg_mapping(self, fname):
        alberty_name_to_cid = {}
        csv_file = csv.reader(open(fname, 'r'))
        csv_file.next()
        for row in csv_file:
            if (row[0] == "" or row[2] == ""):
                continue
            cid = int(row[0])
            alberty_name = row[2]
            alberty_name_to_cid[alberty_name] = cid
        return alberty_name_to_cid
    
    def __init__(self):
        alberty_name_to_pmap = self.read_alberty_mathematics("../data/alberty_mathematica.txt")
        alberty_name_to_cid = self.read_alberty_kegg_mapping("../data/alberty_names.csv")
        self.cid2pmap = {}
        for name in sorted(alberty_name_to_cid.keys()):
            cid = alberty_name_to_cid[name]
            self.cid2pmap[cid] = alberty_name_to_pmap[name]
    
    def get_cid2pmap(self):
        return self.cid2pmap
        
    def cid_to_dG0(self, cid, pH=default_pH, I=default_I, T=default_T):
        try:
            return pmap_to_dG0(self.cid2pmap[cid], pH, I, T, most_abundant=False)
        except KeyError:
            raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)
    
    def reaction_to_dG0(self, sparse_reaction, pH=default_pH, I=default_I, T=default_T):
        dG0 = 0 # calculate the predicted dG0_r
        for (cid, coeff) in sparse_reaction.iteritems():
            dG0_f = self.cid_to_dG0(cid, pH, I, T)
            dG0 += dG0_f * coeff
        return dG0
    
    def display_pmap(self, cid):
        for ((nH, z), dG0) in self.cid2pmap[cid].iteritems():
            print "C%05d | %2d | %6.2f" % (cid, z, dG0)
    
    def write_data_to_csv(self, csv_fname='../res/alberty.csv'):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'nH', 'charge', 'dG0'])
        for cid in sorted(self.cid2pmap.keys()):
            for ((nH, z), dG0) in self.cid2pmap[cid].iteritems():
                writer.writerow([cid, nH, z, dG0])
                
    def write_transformed_data_to_csv(self, csv_fname='../res/alberty_transformed.csv', pH=default_pH, I=default_I, T=default_T):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'pH', 'I', 'T', 'dG0_tag'])
        for cid in sorted(self.cid2pmap.keys()):
            dG0_tag = self.cid_to_dG0(cid, pH, I, T)
            writer.writerow([cid, pH, I, T, dG0_tag])
    
if (__name__ == '__main__'):
    A = Alberty()
    A.write_data_to_csv()
    A.write_transformed_data_to_csv()
    I = 0.25; T = 300;
    sparse_reaction = {13:-1, 1:-1, 9:2}
    sys.stdout.write("The dG0_r of PPi + H20 <=> 2 Pi: \n\n")
    sys.stdout.write("%5s | %5s | %6s | %6s\n" % ("pH", "I", "T", "dG0_r"))
    for pH in arange(5, 9.01, 0.25):
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % (pH, I, T, A.reaction_to_dG0(sparse_reaction, pH, I=I, T=T)))
