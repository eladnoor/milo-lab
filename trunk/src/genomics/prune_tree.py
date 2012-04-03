#!/usr/bin/python

import csv
import dendropy
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
                          help="Annotations to use for pruning")
    return opt_parser


def MakeBool(val):
    return val.lower() == 'true'    


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
    labels_to_prune = set()
    or_lambda = lambda x,y: x or y
    print 'Tree has', len(leaves), 'leaf nodes'    
    for leaf in leaves:
        taxon = leaf.taxon
        taxon_id = taxon.label
        
        d = leaf.comment_metadata
        annot_values = []
        for annot in annot_names:
            value = MakeBool(d.get(annot, 'False'))
            setattr(leaf, annot, value)
            leaf.annotate(annot)
            annot_values.append(value)
        
        val = reduce(or_lambda, annot_values)
        if not val:
            labels_to_prune.add(taxon_id)
    
    print 'Will prune', len(labels_to_prune), 'leaves'
    tree.prune_taxa_with_labels(labels_to_prune)
    tree.write_to_path(options.output_filename, 'nexus')

    bool_colormap = color.ColorMap([True, False])
    for annot in annot_names:
        fname = annot + '.csv'
        f = open(fname, 'w')
        w = csv.writer(f)
        for leaf in leaves:
            taxon = leaf.taxon
            taxon_id = taxon.label
            d = leaf.comment_metadata
            value = MakeBool(d.get(annot, 'False'))
            w.writerow([taxon_id, bool_colormap[value], value])
        
        f.close()
            
            
        
    
    
if __name__ == '__main__':
    Main()