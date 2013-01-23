######################################################
#
# HashedMolecule is a class that defines uniqueness 
# of smiles molecules.
# This class provides comparison (eq) and hashing for
# smiles.
#
# The main script uses this class by reducing the number
# of a given smiles list into the same list without
# repetitions. If given three arguments, the script
# return only the difference molecules from the first
# list (i.e. only the ones that does not appear in the
# second list).
#
######################################################

from openbabel import OBConversion, OBMol


class HashedMolecule (object):
    
    def __init__ (self, smile):
        self._smile = smile.strip()
        conv = OBConversion()
        if not conv.SetInAndOutFormats('smi','inchi'):
            raise 'Problem with openbabel'
        mol = OBMol()
        if not conv.ReadString(mol,self._smile):
            raise TypeError, "No such smile: %s"%self._smile
        self._inchi = conv.WriteString(mol).strip()
        
    def __hash__ (self):
        return self._inchi.__hash__()

    def __eq__ (self, other):
        return self._inchi == other._inchi

    def __repr__ (self):
        return self._smile

    __str__ = __repr__
    
    def getInchi (self):
        return self._inchi




if __name__ == "__main__":
    from sys import argv
    if len(argv) < 3:
        print "Usage: %s <list_file> [<list_difference>] <output>"%argv[0]
        exit(1)
    
    infile = open(argv[1])
    mol_set = set()
    for smile in infile:
        if smile != '\n':
            mol_set.add(HashedMolecule(smile))
    infile.close()
    
    if len(argv) == 4:
        infile = open(argv[2])
        diff_set = set()
        for smile in infile:
            if smile != '\n':
                diff_set.add(HashedMolecule(smile))
        infile.close()
        mol_set = mol_set-diff_set
    
    outfile = open(argv[-1],'w')
    outfile.write('\n'.join(map(lambda x: str(x), mol_set)))
    
