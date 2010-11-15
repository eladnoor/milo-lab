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
G.load_cid2pmap(recalculate=False)

if False:
    H = Hatzi()
    (pH, I, T) = (7,0.1,298.15)
    
    cids = []
    smiles = []

    
    #cids = [65]
    smiles = ["C(O)(=O)CC(=O)O"]
    
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

def calc_pKa(group_0, group_1):
    glist = [G.all_groups[i] for i in G.nonzero_groups]
    index = [glist.index(group_0), glist.index(group_1)]
    dG0_f = [G.group_contributions[i] for i in index]
    return (dG0_f[0] - dG0_f[1])/(R*default_T*log(10))

if True: # calculate the pKa for some common groups
    
    print "-NH2 (0 -> 1)", calc_pKa((u"-N", 2, 0), (u"-N", 3, 1))
    print "-COO (-1 -> 1)", calc_pKa((u"-COO", 0, -1), (u"-COO", 1, 0))
    print "-OPO3 (-2 -> -1)", calc_pKa((u"-OPO3", 0, -2), (u"-OPO3", 1, -1))
    print "-CO-OPO3 (-2 -> -1)", calc_pKa((u"CO-OPO3", 0, -2), (u"CO-OPO3", 1, -1))
    print "-CO-OPO3 (-1 -> 0)", calc_pKa((u"CO-OPO3", 1, -1), (u"CO-OPO3", 2, 0))
