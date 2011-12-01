from pygibbs.kegg import Kegg
from toolbox.html_writer import HtmlWriter
from pygibbs.kegg_errors import KeggParseException
from toolbox.molecule import OpenBabelError

if __name__ == "__main__":
    kegg = Kegg.getInstance()
    
    
    graph = {}
    for rid in kegg.get_all_rids():
        r = kegg.rid2reaction(rid)
        for cid1 in r.sparse.keys():
            for cid2 in r.sparse.keys():
                if r.sparse[cid1] * r.sparse[cid2] < 0:
                    graph.setdefault(cid1, set()).add(cid2)
    
    queue = [355]
    cofactors = set([1,2,3,4,5,6,7,8,9,10,11,13,14,20,28,30])
    html_writer = HtmlWriter('../res/kegg_bfs.html')
    
    for i in xrange(3):
        next_queue = set()
        cofactors.update(queue)
        while queue:
            cid = queue.pop(0)
            next_queue.update(graph[cid])
        queue = list(next_queue.difference(cofactors))
        
        for cid in queue:
            try:
                html_writer.write(kegg.cid2mol(cid).ToSVG())
                html_writer.write(kegg.cid2name(cid))
            except (KeggParseException, OpenBabelError):
                html_writer.write(kegg.cid2name(cid))
            
