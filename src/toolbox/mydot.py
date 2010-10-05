import pydot
import xdot, gtk

def dict2Dot(dict):
    """
        Converts a dictionary to a directed Dot graph.
        In this dictionary, the keys are the nodes and the values are lists of children nodes.
        If a node appears as a neighbor, but is not a key in the dictionary, it is assumed to be a leaf (i.e. has no children).
    """
    Gdot = pydot.Dot()
    
    all_nodes = set()
    for (node, children) in dict.iteritems():
        all_nodes.add(node)
        all_nodes |= set(children)
    
    for node in all_nodes:
        Gdot.add_node(pydot.Node(node, None))
    
    for (node, children) in dict.iteritems():
        for child in children:
            Gdot.add_edge(pydot.Edge(node, child, None))
            
    return Gdot


if (__name__ == "__main__"):
    
    G = {'a' : ['b', 'c', 'd'], 'b' : ['c', 'd'], 'd' : ['a', 'e']}
    Gdot = dict2Dot(G)
    
    win = xdot.DotWindow()
    win.connect('destroy', gtk.main_quit)
    win.set_filter('dot')
    fname = '.dot'
    Gdot.write(fname, format='dot')
    win.open_file(fname)
    gtk.main()