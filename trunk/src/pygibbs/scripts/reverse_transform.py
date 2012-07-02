from toolbox.database import SqliteDatabase
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.nist import NistRowData
import numpy as np
from pygibbs.thermodynamic_constants import R
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg import Kegg
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics

############################################################################

db = SqliteDatabase("../res/gibbs.sqlite")
kegg = Kegg.getInstance()
dissociation = DissociationConstants.FromPublicDB()

nist_row_data = NistRowData()
nist_row_data.pH = 7.77
nist_row_data.I = 0.722 
nist_row_data.T = 298.15
nist_row_data.pMg = 14
nist_row_data.evaluation = "A"
nist_row_data.url = ""
nist_row_data.ref_id = ""
#nist_row_data.reaction = Reaction(["triose isomerase"], {111:-1, 118:1})
#cid2nH = {111:6, 118:5}
#cid2min_dG0 = {111:-1296.3, 118:-1288.6}

nist_row_data.reaction = Reaction(["triose isomerase"], {1:-1, 78:-1, 22:1, 14:1, 463:1})
train_species = PsuedoisomerTableThermodynamics.FromCsvFile(
            '../data/thermodynamics/dG0.csv')
cid2nH = {}
cid2dG0 = {}
for cid in train_species.get_all_cids():
    pmap = train_species.cid2PseudoisomerMap(cid)
    pmatrix = pmap.ToMatrix()
    if len(pmatrix) != 1:
        raise Exception("C%05d has more than one species in the training set" % cid)
    cid2nH[cid] = pmatrix[0][0] # ToMatrix returns tuples of (nH, z, nMg, dG0)
    cid2dG0[cid] = pmatrix[0][3]

ddG0_r = 0
dG0_r = 0
dG0_tag_r = 0
for cid, coeff in nist_row_data.reaction.sparse.iteritems():
    diss = dissociation.GetDissociationTable(cid)
    diss.SetFormationEnergyByNumHydrogens(dG0=cid2dG0[cid], nH=cid2nH[cid], nMg=0)
    for pdata in diss.GenerateAll():
        pdata.name = kegg.cid2name(cid)
        print pdata
        if pdata.hydrogens == cid2nH[cid] and pdata.magnesiums == 0:
            dG0_f = pdata.dG0
            
    dG0_f_tag = diss.Transform(pH=nist_row_data.pH, I=nist_row_data.I, 
                               pMg=nist_row_data.pMg, T=nist_row_data.T)
    
    ddG0_f = diss.GetDeltaDeltaG0(
        pH=nist_row_data.pH, I=nist_row_data.I, pMg=nist_row_data.pMg, 
        T=nist_row_data.T, nH=cid2nH[cid])
    
    print "dG'0_f = %.1f" % dG0_f_tag
    print "dG0_f (obs) = %.1f" %  dG0_f
    print "ddG0_f = %.1f" % ddG0_f
    dG0_r += coeff * dG0_f
    dG0_tag_r += coeff * dG0_f_tag
    ddG0_r += coeff * ddG0_f

# this is actually dG0'_r:
nist_row_data.dG0_r = dG0_tag_r
nist_row_data.K_tag = np.exp(-nist_row_data.dG0_r/(R*nist_row_data.T))

print '-'*50

print "dG0'_r = %.1f" % dG0_tag_r
print "dG0_r  = %.1f" % dG0_r
print "dG0'_r - dG0_r = %.1f" % (dG0_tag_r - dG0_r)
print "ddG0_r = %.1f" % ddG0_r

print '-'*50

data = dissociation.ReverseTransformNistRows([nist_row_data], cid2nH=cid2nH)
print "dG'0_r = %.1f" % data['dG0_r_tag'][0, 0]
print "dG0_r = %.1f" % data['dG0_r'][0, 0]
print "ddG0_r = %.1f" % data['ddG0_r'][0, 0]

print '-'*50
