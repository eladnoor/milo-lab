from toolbox.ods import ODSDictReader
import pybel
from toolbox.smarts_util import FindSmarts
import logging

class EnzymeClass(object):
    
    def __init__(self, ec, substrate, product):
        self.ec = ec
        self.substrate = substrate
        self.product = product
        try:
            self.smarts_subs = pybel.Smarts(substrate)
            self.smarts_prod = pybel.Smarts(product)
        except IOError:
            logging.warning("failed parsing: %s >> %s" % (substrate, product))
            self.smarts_subs = None
            self.smarts_prod = None
        
    def __str__(self):
        return "%9s : %s => %s" % (self.ec, self.substrate, self.product)

    @staticmethod
    def FindSmarts(mol, smarts):
        """
        Corrects the pyBel version of Smarts.findall() which returns results as tuples,
        with 1-based indices even though Molecule.atoms is 0-based.
    
        Args:
            mol: the molecule to search in.
            smarts_str: the SMARTS query to search for.
        
        Returns:
            The re-mapped list of SMARTS matches.
        """
        shift_left = lambda m: [(n - 1) for n in m] 
        return map(shift_left, smarts.findall(mol))    
    
    def React(self, mol):
        if self.smarts_subs and self.smarts_prod:
            products = []
            for atoms in EnzymeClass.FindSmarts(mol, self.smarts_subs):
                products.append(atoms)
            return products
        else:
            return []

class EnzymeMarketplace(object):
    
    def __init__(self):
        self.FromFile()
    
    def FromFile(self, fname='../../data/pathfinder/ReBiT.ods'):
        self.enzyme_class_dict = {}
        for row_dict in ODSDictReader(fname):
            id = int(float(row_dict['EnzymeID']))
            ec = row_dict['EC Number']
            
            substrate = row_dict['Prim Substr Func Grp']
            subs_secondary = row_dict['Sec Substr Func Grp']
            if subs_secondary:
                substrate = '(%s).(%s)' % (substrate, subs_secondary)
            substrate = str(substrate)
            
            product = row_dict['Prim Prod Func Grp']
            prod_secondary = row_dict['Sec Prod Func Grp']
            if prod_secondary:
                product = '(%s).(%s)' % (product, prod_secondary)
            product = str(product)  
    
            self.enzyme_class_dict[id] = EnzymeClass(ec, substrate, product)
        
    def React(self, mol):
        products = []
        for id, eclass in sorted(self.enzyme_class_dict.iteritems()):
            for product in eclass.React(mol):
                products.append((id, product))
        return products

def main():
    mol = pybel.readstring('smiles', 'C(O)(=O)C(=O)O')
    emp = EnzymeMarketplace()
    print emp.React(mol)
    
if __name__ == "__main__":
    main()