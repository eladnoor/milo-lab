######################################
#
# This script runs over all the cids
# and returns their smiles as a list.
# Creates a local file named "smiles.txt"
# which contains the list of smiles.
#
# Note! This script returns segmentation
# fault if running via obable 2.3.2
#
# This script takes no arguments and should
# run from the src directory.
#
######################################

from pygibbs.kegg import Kegg
from openbabel import OBConversion, OBMol
from sys import argv

kegg = Kegg.getInstance()
cids = kegg.get_all_cids()
print 'Found %d cids!'%len(cids)

smiles = open('smiles.txt', 'w')

conv = OBConversion()
if not conv.SetInAndOutFormats('inchi','smi'):
    raise 'Problem with openbabel'

for cid in cids:
    try:
        inchi = kegg.cid2inchi(cid)
        mol = OBMol()
        if not conv.ReadString(mol, inchi):
            raise Exception
        smiles.write(conv.WriteString(mol).strip() + '\n')
    except:
        print "C%05d cannot be converted to smiles" % cid
    
smiles.close()
