import sys
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError, GroupDecompositionError
from hatzimanikatis import Hatzi
from pygibbs.thermodynamics import Thermodynamics, R, default_T
import pybel
from pygibbs.kegg import KeggParseException
from pylab import log
        
G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="pathways")
G.read_compound_abundance("../data/thermodynamics/compound_abundance.csv")
G.init()
G.load_cid2pmap()
G.write_data_to_json("../res/group_contribution.json", G.kegg())

if False:
    H = Hatzi()
    (pH, I, T) = (7,0.1,298.15)
    
    cids = []
    smiles = []

    
    #cids = [65]
    smiles = ["OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C(O)C(O)C(O)C(O)C(O)C(=O)"]
    
    mols = []
    for cid in cids:
        try:
            mols += [G.kegg().cid2mol(cid)]
        except KeggParseException:
            continue
        except KeyError:
            continue

    for s in smiles:
        mol = pybel.readstring('smiles', s)
        mol.removeh()
        mols += [mol]

    #print pybel.Smarts("*=[N,n;H0;D3;R1;+]").findall(mols[0])
    #mols[0].draw()
    #print pybel.Smarts("[C;H1;R1]").findall(mols[0])
    #sys.exit(0)
        
    for m in mols:
        try:
            m.removeh()
            print G.analyze_decomposition(m)
            
            pmap_m = G.estimate_pmap(m)
            #pmap_h = H.cid2pmap(cid)
            
            print Thermodynamics.pmap_to_table(pmap_m)
            print "dG0'_f (M): ", Thermodynamics.pmap_to_dG0(pmap_m, pH, I, T)
            #print "dG0'_f (H): ", Thermodynamics.pmap_to_dG0(pmap_h, pH, I, T)
        except KeyError as e:
            sys.stderr.write(str(e) + "\n")
        except GroupMissingTrainDataError as e:
            sys.stderr.write(str(e) + "\n")
            for miss_list in e.missing_groups:
                sys.stderr.write(str(miss_list) + "\n")
        except GroupDecompositionError as e:
            sys.stderr.write(str(e) + "\n")
        m.draw()

if False:
    G.analyze_pathway("../data/thermodynamics/pathways.txt")

def calc_pKa_group(group_0, group_1):
    glist = [G.all_groups[i] for i in G.nonzero_groups]
    index = [glist.index(group_0), glist.index(group_1)]
    dG0_f = [G.group_contributions[i] for i in index]
    pKa = (dG0_f[0] - dG0_f[1])/(R*default_T*log(10))
    return "%s, %d (%.1f) -> %d (%.1f), pKa = %.1f" % (group_0[0], group_0[2], dG0_f[0], group_1[2], dG0_f[1], pKa)

def calc_pKa_compound(cid, pseudo_0, pseudo_1):
    pmap = G.cid2pmap_obs[cid]
    dG0_f = [pmap[pseudo_0][0], pmap[pseudo_1][0]]
    return (dG0_f[0] - dG0_f[1])/(R*default_T*log(10))

if True: # calculate the pKa for some common groups
    
    print calc_pKa_group((u"-N", 2, 0), (u"-N", 3, 1))
    print calc_pKa_group((u"-COO", 0, -1), (u"-COO", 1, 0))
    print calc_pKa_group((u"-OPO3", 0, -2), (u"-OPO3", 1, -1))
    print calc_pKa_group((u"CO-OPO3", 0, -2), (u"CO-OPO3", 1, -1))
    print calc_pKa_group((u"CO-OPO3", 1, -1), (u"CO-OPO3", 2, 0))

    #print "%30s: %5.2f" % ("acetyl-P (-2 -> -1)", calc_pKa_compound(227, (3, -2), (4, -1)))
    #print "%30s: %5.2f" % ("acetyl-P (-1 -> 0)", calc_pKa_compound(227, (4, -1), (5, 0)))
    #print "%30s: %5.2f" % ("carbamoyl-P (-2 -> -1)", calc_pKa_compound(169, (2, -2), (3, -1)))
    #print "%30s: %5.2f" % ("carbamoyl-P (-1 -> 0)", calc_pKa_compound(169, (3, -1), (4, 0)))

    #print "%30s: %5.2f" % ("ATP (-4 -> -3)", calc_pKa_compound(2, (12, -4), (13, -3)))
    #print "%30s: %5.2f" % ("ATP (-3 -> -2)", calc_pKa_compound(2, (13, -3), (14, -2)))
    #print "%30s: %5.2f" % ("BPG (-4 -> -3)", calc_pKa_compound(236, (4, -4), (5, -3)))
