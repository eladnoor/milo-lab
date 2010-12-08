import csv
from pylab import log, sqrt, array
from toolbox.util import log_sum_exp
import pseudoisomer

R = 8.31e-3 # kJ/(K*mol)
J_per_cal = 4.184
default_T = 298.15 # K
default_I = 0.1 # mM
default_pH = 7.0
default_c0 = 1 # M
default_pMg = 3


class MissingCompoundFormationEnergy(Exception):
    def __init__(self, value, cid=0):
        self.value = value
        self.cid = cid
    def __str__(self):
        return repr(self.value)    

class Thermodynamics:
    def __init__(self):
        self.pH = default_pH
        self.I = default_I
        self.T = default_T
        self.pMg = default_pMg
        
        self.c_mid = 1e-3
        self.c_range = (1e-6, 1e-2)
        self.bounds = {}
        self.source_string = "Unknown"

    @staticmethod
    def debye_huckel(I):
        return (2.91482 * sqrt(I)) / (1 + 1.6 * sqrt(I))

    @staticmethod
    def correction_function(nH, nMg, z, pH, pMg, I, T):
        """
            nH and z - are the species parameters (can be vectors)
            pH and I - are the conditions, must be scalars
            returns the correction element used in the transform function
            
        Returns:
            The correction, in units of RT.
        """
        DH = Thermodynamics.debye_huckel(I) / (R*T)
        return -nMg * (log(10)*pMg) - nH * (log(10)*pH + DH) + (z**2) * DH

    @staticmethod
    def transform(dG0, nH, z, pH, I, T):
        return dG0 - R*T*Thermodynamics.correction_function(nH, z, pH, I, T)

    @staticmethod
    def array_transform(dG0, nH, nMg, z, pH, pMg, I, T):
        """
            dG0, nH and z - are the species parameters (can be vectors)
            pH and I - are the conditions, must be scalars
            returns the transformed gibbs energy: dG0'
        """
        return (-(R*T) * log_sum_exp(dG0 / (-R*T) +
                Thermodynamics.correction_function(nH, nMg, z, pH, pMg, I, T)))

    def cid2pmap(self, cid):
        raise Exception("method not implemented")

    def cid2pmatrix(self, cid):
        return self.cid2pmap(cid).ToMatrix()

    def get_all_cids(self):
        raise Exception("method not implemented")
        
    def cid_to_dG0(self, cid, pH=None, pMg=None, I=None, T=None):
        pH = pH or self.pH
        I = I or self.I
        T = T or self.T
        pMg = pMg or self.pMg
        return self.cid2pmap(cid).Transform(pH, pMg, I, T, most_abundant=False)
    
    def reaction_to_dG0(self, sparse_reaction, pH=None, I=None, T=None):
        """
            calculate the predicted dG0_r
        """
        return sum([coeff * self.cid_to_dG0(cid, pH, I, T) for (cid, coeff) in sparse_reaction.iteritems()])
    
    def cid_to_bounds(self, cid, use_default=True):
        (curr_c_min, curr_c_max) = self.bounds.get(cid, (None, None))
        if (curr_c_min == None and use_default):
            curr_c_min = self.c_range[0]
        if (curr_c_max == None and use_default):
            curr_c_max = self.c_range[1]
        return (curr_c_min, curr_c_max)

    @staticmethod
    def pmap_to_table(pmap, pH=default_pH, I=default_I, T=default_T):
        s = ""
        s += "%2s | %2s | %7s | %7s\n" % ("nH", "z", "dG0_f", "dG0'_f")
        s += "-" * 35 + "\n"
        for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(pmap):
            s += "%2d | %2d | %7.1f | %7.1f\n" % (nH, z, dG0, Thermodynamics.transform(dG0, nH, z, pH, I, T))
        return s     

    def display_pmap(self, cid):
        for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(self.cid2pmap(cid)):
            print "C%05d | %2d | %2d | %6.2f" % (cid, nH, z, dG0)
    
    def write_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'nH', 'charge', 'nMg', 'dG0'])
        for cid in self.get_all_cids():
            for (nH, z, nMg, dG0) in self.cid2pmap(cid).ToMatrix():
                writer.writerow([cid, nH, z, nMg, dG0])

    def write_data_to_json(self, json_fname, kegg):
        import json

        formations = []
        for cid in self.get_all_cids():
            h = {}
            h["inchi"] = kegg.cid2inchi(cid)
            h["source"] = self.source_string
            h["species"] = []
            for (nH, z, nMg, dG0) in self.cid2pmap(cid).ToMatrix():
                h["species"].append({"nH":nH, "z":z, "nMg":nMg, "dG0_f":dG0})
            formations.append(h)

        json_file = open(json_fname, 'w')
        json_file.write(json.dumps(formations))
        json_file.close()
                
    def write_transformed_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['CID', 'pH', 'pMg', 'I', 'T', 'dG0_tag'])
        for cid in self.get_all_cids():
            dG0_tag = self.cid_to_dG0(cid)
            writer.writerow([cid, self.pH, self.pMg, self.I, self.T, dG0_tag])
            
    def save_energies_to_db(self, db, table_name):
        db.CreateTable(table_name, "cid INT, dG0_f REAL, nH INT, z INT, nMg INT, anchor BOOL")
        for cid in self.get_all_cids():
            for (nH, z, nMg, dG0) in self.cid2pmap(cid).ToMatrix():
                db.Insert(table_name, [cid, dG0, nH, z, nMg, cid in self.anchors])
        db.Commit()

    def load_energies(self, db, table_name):
        self.cid2pmap_dict = {}
        self.anchors = set()
        for row in db.execute("SELECT * FROM %s" % table_name):
            (cid, dG0, nH, z, nMg, anchor) = row
            self.cid2pmap_dict.setdefault(cid, pseudoisomer.PseudoisomerMap())
            self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)
            if (anchor):
                self.anchors.add(cid)
        self.update_cache(self)