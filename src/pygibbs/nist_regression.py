import logging, sys
from pygibbs.nist import Nist
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from toolbox.util import _mkdir
from toolbox.html_writer import HtmlWriter
from pygibbs.group_decomposition import GroupDecomposer

if (__name__ == "__main__"):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    _mkdir("../res/nist")
    
    kegg = Kegg()
    nist = Nist(kegg)
    html_writer = HtmlWriter("../res/nist/regression.html")
    db = SqliteDatabase('../res/gibbs.sqlite')
    
    group_decomposer = GroupDecomposer.FromDatabase(db)
    
    dissociation = DissociationConstants(db, html_writer, kegg)
    dissociation.LoadValuesToDB('../data/thermodynamics/pKa_with_cids.csv')
    cid2pKa_list = dissociation.GetAllpKas()
    
    cids_with_pKa = set(cid2pKa_list.keys())
    cids_in_nist = set(nist.cid2count.keys())
    
    print "CIDs with pKa: ", len(cids_with_pKa)
    print "CIDs in NIST: ", len(cids_in_nist)
    print "CIDs in NIST with pKas: ", len(cids_in_nist.intersection(cids_with_pKa))
    
    print "CIDs in NIST wihtout pKas: "
    for cid in sorted(cids_in_nist):
        try:
            mol = kegg.cid2mol(cid)
            decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
        except Exception:
            print "C%05d - cannot decompose" % cid
            continue
        
        if len(decomposition.PseudoisomerVectors()) > 1:
            print "C%05d - should have pKas" % cid
        else:
            print "C%05d - doesn't have pKas" % cid
            
    for cid in sorted(cids_in_nist.difference(cids_with_pKa)):
        print "C%05d" % cid
    
    html_writer.close()