import sys
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError, GroupDecompositionError
from hatzimanikatis import Hatzi
from pygibbs.thermodynamics import Thermodynamics
import pybel
from pygibbs.kegg import KeggParseException
        
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


if True: # calculate the pKa for some common groups
    glist = [(group_name, charge) for (gid, group_name, protons, charge, smarts_str, focal_atoms) in G.list_of_groups]
    
    i0 = G.nonzero_groups[glist.index(("-COO", 0))]
    i1 = G.nonzero_groups[glist.index(("-COO", -1))]
    
    print G.group_contributions[i0], G.group_contributions[i1]
