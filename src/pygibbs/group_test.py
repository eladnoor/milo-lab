import sys
from groups import GroupContribution, GroupMissingTrainDataError, GroupDecompositionError
        
G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_test")
G.read_compound_abundance("../data/compound_abundance.csv")
G.write_gc_tables()
G.init()

if False:
    for cid in [101]:
        m = G.kegg().cid2mol(cid)
        #print Smarts('[NH1]=[C,c,N,n]').findall(m)
        try:
            print G.analyze_decomposition_cid(cid)
            print G.estimate_dG0_keggcid(cid)
        except KeyError as e:
            sys.stderr.write(str(e) + "\n")
        except GroupMissingTrainDataError as e:
            sys.stderr.write(str(e) + "\n")
            for miss_list in e.missing_groups:
                sys.stderr.write(str(miss_list) + "\n")
        except GroupDecompositionError as e:
            sys.stderr.write(str(e) + "\n")
        m.draw()

G.analyze_pathway("../data/pathways.txt", "../res/pathways.html")
