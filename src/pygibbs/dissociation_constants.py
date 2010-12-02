import csv, re, sys
from pylab import arange, NaN, isnan
from kegg import Kegg

class DissociationConstants:
    def __init__(self):
        self.kegg = Kegg()
    
    def read_csv(self, fname):
        csv_reader = csv.reader(open(fname, 'r'))
        csv_reader.next() # skip title row
        
        last_formula = None
        last_name = None
        for formula, name, step, T, pKa in csv_reader:
            if (name != last_name or formula != last_formula):
                cid = self.find_cid(formula, name)
                
            formula = formula or last_formula
            name = name or last_name
            step = int(step or 1)
            
            print cid, formula, name, step, T, pKa            
            
            last_formula = formula
            last_name = name
    
    def formula_to_atom_bag(self, formula):
        if (formula == None or formula.find("(") != -1 or formula.find(")") != -1 or formula.find("R") != -1):
            raise Exception("non-specific compound formula: " + formula)
        
        atom_bag = {}
        for (atom, count) in re.findall("([A-Z][a-z]*)([0-9]*)", formula):
            if (count == ''):
                count = 1
            else:
                count = int(count)
            atom_bag[atom] = count

        return atom_bag
    
    def find_cid(self, formula, name):
        cid = self.kegg.name2cid(name.lower(), fuzziness_cutoff=1)
        if cid:
            print ">",
        else:
            cid = self.kegg.name2cid(name.lower(), fuzziness_cutoff=0.5)
        if not cid:
            return None
        atom_bag = self.formula_to_atom_bag(formula)
        kegg_atom_bag = self.kegg.cid2atom_bag(cid)
        if (atom_bag != kegg_atom_bag):
            return None
        else:
            return cid
    
if (__name__ == '__main__'):
    dc = DissociationConstants()
    dc.read_csv('../data/thermodynamics/pKa.csv')