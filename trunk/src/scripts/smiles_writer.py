from pygibbs.kegg import Kegg
from openbabel import OBConversion, OBMol
from sys import argv

kegg = Kegg.getInstance()
cids = kegg.get_all_cids()
print 'Found %d cids!'%len(cids)

smiles = open('smiles.txt','a')
nums = open('nums.txt','a')

if len(argv) < 2:
    i = 0
else:
    i = int(argv[1])

conv = OBConversion()
if not conv.SetInAndOutFormats('inchi','smi'):
    raise 'Problem with openbabel'

for j,cid in enumerate(cids[i:]):
    try:
        nums.write(str(j+i) + '\n')
        inchi = kegg.cid2inchi(cid)
        mol = OBMol()
        if not conv.ReadString(mol,inchi):
            print cid
        smiles.write(conv.WriteString(mol))
        #smiles.write(kegg.cid2smiles(cid) + '\n')
    except:
        print cid
    
smiles.close()
nums.close()
