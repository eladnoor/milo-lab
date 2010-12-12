import sys
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError, GroupDecompositionError
from hatzimanikatis import Hatzi
from pygibbs.thermodynamics import Thermodynamics, R, default_T
import pybel
from pygibbs.kegg import KeggParseException
from pylab import log
from toolbox import database
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import Group
        
db = database.SqliteDatabase('../res/gibbs.sqlite')
html_writer = HtmlWriter('../res/dG0_test.html')
G = GroupContribution(db, html_writer)
G.read_compound_abundance("../data/thermodynamics/compound_abundance.csv")
G.init()
G.load_cid2pmap()
G.write_data_to_json("../res/group_contribution.json", G.kegg())

if False:
    H = Hatzi()
    (pH, pMg, I, T) = (7, 3, 0.1, 298.15)
    
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
            
            print pmap_m
            print "dG0'_f (M): ", pmap_m.Transform(pH, pMg, I, T)
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
    index_full = [G.groups_data.Index(group_0), G.groups_data.Index(group_1)]
    index_short = [G.nonzero_groups.index(i) for i in index_full]
    dG0_f = [G.group_contributions[i] for i in index_short]
    pKa = (dG0_f[0] - dG0_f[1])/(R*default_T*log(10))
    return pKa

def calc_pKa_compound(cid, pseudo_0, pseudo_1):
    pmap = G.cid2pmap_obs[cid]
    dG0_f = [pmap[pseudo_0][0], pmap[pseudo_1][0]]
    return (dG0_f[0] - dG0_f[1])/(R*default_T*log(10))

if True: # calculate the pKa for some common groups
    
    print "-NH2 -> -NH3[+] : %.1f" % calc_pKa_group(Group(None, "-N", 2, 0, 0), Group(None, "-N", 3, 1, 0))
    print "-COO[-] -> -COOH : %.1f" % calc_pKa_group(Group(None, "-COO", 0, -1, 0), Group(None, "-COO", 1, 0, 0))
    print "-OPO3[2-] -> -HOPO3[-] : %.1f" % calc_pKa_group(Group(None, "-OPO3", 0, -2, 0), Group(None, "-OPO3", 1, -1, 0))
    print "-CO-OPO3[2-] -> -CO-HOPO3[-] : %.1f" % calc_pKa_group(Group(None, "CO-OPO3", 0, -2, 0), Group(None, "CO-OPO3", 1, -1, 0))
    print "-CO-HOPO3[-] -> -CO-H2OPO3 : %.1f" % calc_pKa_group(Group(None, "CO-OPO3", 1, -1, 0), Group(None, "CO-OPO3", 2, 0, 0))

    print "-OPO3.Mg -> -OPO3[2-] : %.1f" % calc_pKa_group(Group(None, "-OPO3", 0, -2, 0), Group(None, "-OPO3", 0, 0, 1))

    #print "%30s: %5.2f" % ("acetyl-P (-2 -> -1)", calc_pKa_compound(227, (3, -2), (4, -1)))
    #print "%30s: %5.2f" % ("acetyl-P (-1 -> 0)", calc_pKa_compound(227, (4, -1), (5, 0)))
    #print "%30s: %5.2f" % ("carbamoyl-P (-2 -> -1)", calc_pKa_compound(169, (2, -2), (3, -1)))
    #print "%30s: %5.2f" % ("carbamoyl-P (-1 -> 0)", calc_pKa_compound(169, (3, -1), (4, 0)))

    #print "%30s: %5.2f" % ("ATP (-4 -> -3)", calc_pKa_compound(2, (12, -4), (13, -3)))
    #print "%30s: %5.2f" % ("ATP (-3 -> -2)", calc_pKa_compound(2, (13, -3), (14, -2)))
    #print "%30s: %5.2f" % ("BPG (-4 -> -3)", calc_pKa_compound(236, (4, -4), (5, -3)))
