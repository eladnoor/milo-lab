import sys
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError
from pygibbs.group_decomposition import GroupDecompositionError
import pybel
from pygibbs.kegg import KeggParseException
from toolbox import database
from toolbox.html_writer import HtmlWriter

def main(G):        
    G.init()
    G.load_cid2pmap()
    
    (pH, pMg, I, T) = (7, 3, 0.1, 298.15)
    
    cids = []
    smiles = ["C(=O)(O)C(=O)CC(=O)O"]
    smarts = pybel.Smarts("[C;H1,H2][C;H0](=O)C")
    
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
    
    for m in mols:
        try:
            m.removeh()
            print smarts.findall(mol)
            print G.analyze_decomposition(m)
            pmap_m = G.estimate_pmap(m)
            print pmap_m
            print "dG0'_f (M): ", pmap_m.Transform(pH, pMg, I, T)
        except KeyError as e:
            sys.stderr.write(str(e) + "\n")
        except GroupMissingTrainDataError as e:
            sys.stderr.write(str(e) + "\n")
            for miss_list in e.missing_groups:
                sys.stderr.write(str(miss_list) + "\n")
        except GroupDecompositionError as e:
            #sys.stderr.write(str(e) + "\n")
            pass
        m.draw()
    
    if True: # calculate the pKa for some common groups
        print "-NH2 -> -NH3[+] : pKa = %.1f" % G.calc_pKa_group("-N", 3, 1)
        print "-COO[-] -> -COOH : pKa = %.1f" % G.calc_pKa_group("-COO", 1, 0)
        print "-OPO3[2-] -> -HOPO3[-] : pKa = %.1f" % G.calc_pKa_group("-OPO3", 1, -1)
        print "-CO-OPO3[2-] -> -CO-HOPO3[-] : pKa = %.1f" % G.calc_pKa_group("CO-OPO3", 1, -1, 0)
        print "-CO-HOPO3[-] -> -CO-H2OPO3 : pKa = %.1f" % G.calc_pKa_group("CO-OPO3", 2, 0, 0)
        print "-OPO3.Mg -> -OPO3[2-] : pK_Mg = %.1f" % G.calc_pK_Mg_group("-OPO3", 0, 0, 1)

if __name__ == "__main__":
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/dG0_test.html')
    G = GroupContribution(db, html_writer)
    main(G)
    #G.analyze_pathway("../data/thermodynamics/pathways.txt")
    
