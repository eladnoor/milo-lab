#!/usr/bin/python

import csv
import dendropy
import re
import sys

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
    opt_parser.add_option("-n", "--name_col",
                          dest="name_col",
                          default="Genome Name",
                          help="Name of the column containing the full name ID.")
    opt_parser.add_option("-i", "--taxon_id_col",
                          dest="taxon_id_col",
                          default="NCBI Taxon ID",
                          help="Name of the column containing the taxonomy ID.")
    return opt_parser


def NormalizeName(name):
    tmp = re.sub('_str_', '_', name.lower().strip())
    tmp = re.sub('_subsp_', '_', tmp)
    tmp = tmp.replace('sp.', 'sp')
    return re.sub('[\s-]', '_', tmp)


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.tree_filename and options.genomes_filename
    assert options.output_filename
    
    print 'Reading tree from', options.tree_filename
    tree = dendropy.Tree()
    tree.read_from_path(options.tree_filename, options.tree_format)
        
    print 'Reading genome data from', options.genomes_filename
    f = open(options.genomes_filename)
    r = csv.DictReader(f)
    taxons_by_name = {}
    for row in r:
        name = row.get(options.name_col, None)
        taxon_id = row.get(options.taxon_id_col, None)
        if not taxon_id or not name:
            continue
        
        norm_name = NormalizeName(name)
        if norm_name in taxons_by_name: 
            print 'Duplicate name', name, '(normed: %s)' % (norm_name)
            continue
        
        taxons_by_name.setdefault(norm_name, set()).add(taxon_id)
        
    f.close()
    
    print 'Found', len(taxons_by_name), 'names'

    leaves = tree.leaf_nodes()
    print 'Tree has', len(leaves), 'leaf nodes'    
    unknown_names = 0
    for leaf in leaves:
        taxon = leaf.taxon
        label = taxon.label
        norm_name = NormalizeName(label)
        if norm_name not in taxons_by_name:
            print 'Unknown name', label, '(normed: %s)' % (norm_name)
            unknown_names += 1
    
    print unknown_names, 'unknown names'
        

    #tree.write_to_path(options.output_filename, 'nexus')

    
    
if __name__ == '__main__':
    Main()