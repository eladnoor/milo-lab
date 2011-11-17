#!/usr/bin/python

import csv
import dendropy
import math
import re
import sys

import pylab
import scipy
from matplotlib.font_manager import FontProperties

from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-t", "--tree_filename",
                          dest="tree_filename",
                          help="Phylogenetic tree data")
    opt_parser.add_option("-f", "--tree_format",
                          dest="tree_format",
                          default="newick",
                          help="Newick format")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="output annotated and pruned tree")
    opt_parser.add_option("-g", "--genomes_filename",
                          dest="genomes_filename",
                          help="filename of CSV to use for labeling")
    opt_parser.add_option("-i", "--taxon_id_col",
                          dest="taxon_id_col",
                          default="NCBI Taxon ID",
                          help="Name of the column containing the taxonomy ID.")
    opt_parser.add_option("-c", "--cols",
                          dest="cols",
                          help="Columns to use for labeling")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.tree_filename and options.genomes_filename
    assert options.output_filename
    assert options.cols
    
    print 'Reading tree from', options.tree_filename
    tree = dendropy.Tree()
    tree.read_from_path(options.tree_filename, options.tree_format)
    
    cols = [s.strip() for s in options.cols.split(',')]
    print 'Will annotate leaves with labels:', ', '.join(cols)
    
    print 'Reading genome data from', options.genomes_filename
    f = open(options.genomes_filename)
    r = csv.DictReader(f)
    rows_by_taxon = {}
    for row in r:
        taxon_id = row.get(options.taxon_id_col, None)
        if not taxon_id:
            continue
        
        if taxon_id in rows_by_taxon: 
            print 'Duplicate taxon ID', taxon_id
        
        rows_by_taxon[taxon_id] = row
    f.close()
    
    all_tree_labels = set(tree.taxon_set.labels())
    all_genome_labels = set([t for t in rows_by_taxon.iterkeys()])
    labels_to_prune = all_tree_labels.difference(all_genome_labels)
    print 'Pruning', len(labels_to_prune), 'of', len(all_tree_labels), 'tree nodes'
    tree.prune_taxa_with_labels(labels_to_prune)


    leaves = tree.leaf_nodes()
    print 'Tree has', len(leaves), 'leaf nodes'    
    print 'Annotating leaves with data from columns %s' % ','.join(cols)
    for leaf in leaves:
        taxon = leaf.taxon
        taxon_id = taxon.label
        if taxon_id not in rows_by_taxon:
            print 'Couldn\'t find taxon ID', taxon_id
            continue
        
        row = rows_by_taxon[taxon_id]
        for col in cols:
            val = row.get(col, None)
            attr_name = re.sub('[\s]', '_', col)
            setattr(leaf, attr_name, val)
            leaf.annotate(attr_name)

    tree.write_to_path(options.output_filename, 'nexus')

    
    
if __name__ == '__main__':
    Main()