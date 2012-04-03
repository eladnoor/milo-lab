#!/usr/bin/python

import csv
import dendropy
import logging
import sys

import pylab
from matplotlib.font_manager import FontProperties
from toolbox import color

from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-t", "--tree_filename",
                          dest="tree_filename",
                          help="Phylogenetic tree in newick format")
    opt_parser.add_option("-f", "--tree_format",
                          dest="tree_format",
                          default="nexus",
                          help="Tree file format")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="output annotated and pruned tree")
    opt_parser.add_option("-a", "--annotations",
                          dest="annotations",
                          help="Annotations to make ITOL files for.")
    opt_parser.add_option("-e", "--edge_length_threshold",
                          type="float",
                          default=0.0,
                          dest="edge_length_threshold",
                          help="Edge length threshold.")
    return opt_parser


def MakeBool(val):
    return val.lower() == 'true'    


def MakeBoolDict(d):
    return dict((k, MakeBool(v)) for k,v in d.iteritems())


def OrDicts(d1, d2):
    keyset = set(d1.keys()).union(d2.keys())
    d = {}
    for k in keyset:
        d[k] = d1.get(k, False) or d2.get(k, False)
    return d

def AreLeaves(nodes):
    for n in nodes:
        if not n.is_leaf():
            return False
    return True

def MaybeMergeChildren(parent_node):
    children = parent_node.child_nodes()
    assert len(children) == 2
    if not AreLeaves(children):
        logging.debug('Not both children are leaves. Bailing.')
        return False
    
    # Make the new dictionaries and edge lengths
    child_dicts = [c.d for c in children]
    child_lengths = [c.edge.length for c in children]
    virtual_count = sum(c.count for c in children)
    max_length_idx = pylab.argmax(child_lengths)
    merged_d = reduce(OrDicts, child_dicts)
    label = children[max_length_idx].taxon.label
    
    logging.debug('Merging 2 children with edge lengths %s',
                  child_lengths)
    
    # Remove children and update the parent
    map(parent_node.remove_child, children)
    parent_node.edge.length += child_lengths[max_length_idx]
    parent_node.d = merged_d
    parent_node.count = virtual_count
    for k, v in parent_node.d.iteritems():
        setattr(parent_node, k, v)
        parent_node.annotate(k)
    
    # Set up a taxon for the parent according to the
    # most distinct child.
    # TODO(flamholz): indicate somehow that this was merged.
    taxon = dendropy.Taxon()
    taxon.label = label
    parent_node.taxon = taxon
    
    
    return True
    


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.tree_filename and options.tree_format
    assert options.output_filename
    assert options.annotations
    
    print 'Reading tree from', options.tree_filename
    tree = dendropy.Tree()
    tree.read_from_path(options.tree_filename, options.tree_format,
                        extract_comment_metadata=True)
    
    annot_names = [s.strip() for s in options.annotations.split(',')]
    print 'Will annotate leaves with labels:', ', '.join(annot_names)
    
    leaves = tree.leaf_nodes()
    print 'Tree has %d leaf nodes' % len(leaves)
    
    print 'Pruning leaves'
    for leaf in leaves:
        leaf.d = MakeBoolDict(leaf.comment_metadata)
        leaf.count = 1
        for k, v in leaf.d.iteritems():
            setattr(leaf, k, v)
            leaf.annotate(k)

    while True:
        changed = False
        for e in tree.leaf_edge_iter():
            if (e.tail_node is not None and
                e.length < options.edge_length_threshold):
                changed |= MaybeMergeChildren(e.tail_node)
        
        
        if not changed:
            break
            
        
    tree.write_to_path(options.output_filename, 'nexus',
                       suppress_annotations=True)

    leaves = tree.leaf_nodes()
    print 'Tree now has %d leaf nodes' % len(leaves)
    labels = [False] + ['%s_%s' % (True, a) for a in annot_names]
    bool_colormap = color.ColorMap(labels)
    
    for annot in annot_names:
        fname = annot + '.csv'
        f = open(fname, 'w')
        w = csv.writer(f)
        for leaf in leaves:
            taxon = leaf.taxon
            taxon_id = taxon.label
            d = leaf.annotations()
            value = d.get(annot, (False, None))[0]
            if value is True:
                value = '%s_%s' % (value, annot)
            w.writerow([taxon_id, bool_colormap[value], value])
        
        f.close()
    

    fname = 'Node_Counts.csv'
    f = open(fname, 'w')
    w = csv.writer(f)
    for leaf in leaves:
        taxon = leaf.taxon
        taxon_id = taxon.label
        w.writerow([taxon_id, leaf.count])
    
    f.close()
            
            
        
    
    
if __name__ == '__main__':
    Main()