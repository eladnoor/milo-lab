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
                          help="Phylogenetic tree in newick format")
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


def ParseTaxonID(label):
    return label
    """
    m = re.match('taxon__(\d+)', label)
    g = m.groups()
    if not g:
        return None
    return g[0]
    """

def AllChildrenAreLeaves(node):
    children = node.child_nodes()
    for node in children:
        if not node.is_leaf():
            return False
    return True


def MergeDists(dist, to_add):
    for key in to_add:
        dist[key] = dist.get(key, 0) + to_add[key]


def MaybeCollapseChildren(node):
    if not AllChildrenAreLeaves(node):
        return False
    
    children = node.child_nodes()
    glyc_types = set()
    total_children = 0
    merged_oxygen_req_dist = {}
    for c in children:
        glyc_types.add(c.glyc_type)
        total_children += c.total_children
        MergeDists(merged_oxygen_req_dist, c.oxygen_req_dist)
            
    ret = False
    if len(glyc_types) == 1:
        c0taxon = children[0].taxon
        for child in children:
            node.remove_child(child)
            
        node.glyc_type = glyc_types.pop()
        node.oxygen_req_dist = merged_oxygen_req_dist
        node.total_children = total_children
        taxon = dendropy.Taxon()
        taxon.label = c0taxon.label
        node.taxon = taxon
        node.annotate('glyc_type')
        ret = True
    return ret
    

def MakeBool(val):
    l = val.lower()
    return l == 'true'


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.tree_filename and options.genomes_filename
    assert options.output_filename
    assert options.cols
    
    print 'Reading tree from', options.tree_filename
    tree = dendropy.Tree()
    tree.read_from_path(options.tree_filename, 'newick')
    
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
    
    glyc_types = {}
    for tax_id, row in rows_by_taxon.iteritems():
        all_vals = {}
        for col in cols:
            val = row.get(col, None)
            if not val:
                print 'Failed to retrieve col %s from row %s' % (col, row)
                continue
                        
            all_vals[col] = MakeBool(val)
        
        glyc_type = 'Neither'
        both = reduce(lambda x,y: x and y, all_vals.values())
        if both:
            glyc_type = 'Both'
        else:
            l = [col for col,val in all_vals.iteritems()
                 if val == True]
            if l:
                glyc_type = l[0]
        
        glyc_types[tax_id] = glyc_type
    
    all_tree_labels = set(tree.taxon_set.labels())
    all_genome_labels = set([t for t in rows_by_taxon.iterkeys()])
    labels_to_prune = all_tree_labels.difference(all_genome_labels)
    print 'Pruning', len(labels_to_prune), 'of', len(all_tree_labels), 'tree nodes'
    tree.prune_taxa_with_labels(labels_to_prune)

    neither_type = set([t for t,glyc in glyc_types.iteritems()
                        if glyc == 'Neither'])
    print 'Pruning', len(neither_type), 'more tree nodes'
    tree.prune_taxa_with_labels(neither_type)

    oxygen_reqs = {}
    for tax_id, row in rows_by_taxon.iteritems():
        oxygen_req = row.get('Broad Oxygen Requirement', 'None given')
        oxygen_req = oxygen_req.strip() or 'None given'
        oxygen_reqs[tax_id] = oxygen_req

    leaves = tree.leaf_nodes()
    print 'Tree has', len(leaves), 'leaf nodes'    
    for leaf in leaves:
        taxon = leaf.taxon
        label = taxon.label
        taxon_id = ParseTaxonID(label)
        if taxon_id not in glyc_types:
            print 'Couldn\'t find taxon ID', taxon_id
            continue
        
        oxygen_req = oxygen_reqs[taxon_id]
        leaf.oxygen_req_dist = {'Aerobe': 0,
                                'Anaerobe': 0,
                                'Facultative': 0,
                                'None given': 0}
        leaf.oxygen_req_dist[oxygen_req] = 1
        leaf.glyc_type = glyc_types[taxon_id]
        leaf.total_children = 1
        leaf.annotate('glyc_type')
    
    print 'Collapsing uniform clades'
    while True:
        leaf_parents = set()
        for l in tree.leaf_nodes():
            leaf_parents.add(l.parent_node)
        
        changed = False
        for node in leaf_parents:
            changed |= MaybeCollapseChildren(node)
            
        if not changed:
            break

    leaves = tree.leaf_nodes()
    print 'Tree now has', len(leaves), 'leaf nodes'

    tree.write_to_path(options.output_filename, 'nexus')

    colors = {'None given': '#000000',
              'Neither': '#000000',
              'Anaerobe': '#ff0000',
              cols[0]: '#CC0000',
              'Aerobe': '#00ff00',
              cols[1]: '#009900',
              'Facultative': '#0000ff',
              'Both': '#0033ff'}

    f = open('oxygen_req_dists.csv', 'w')
    order = ['Aerobe', 'Anaerobe', 'Facultative', 'None given']
    w = csv.writer(f)
    w.writerow(['LABELS'] + order)
    w.writerow(['COLORS'] + [colors[i] for i in order])
    for leaf in leaves:
        id = leaf.taxon.label
        d = leaf.oxygen_req_dist
        counts = [d[i] for i in order]
        total = float(sum(counts))
        pcts = [100*float(i)/total for i in counts]
        row = [id] + pcts
        w.writerow(row)
    f.close()
    
    f = open('glyc_types.csv', 'w')
    w = csv.writer(f)
    for leaf in leaves:
        id = leaf.taxon.label
        glyc_type = leaf.glyc_type
        w.writerow([id, colors[glyc_type], glyc_type])
    f.close()
    
    matrix = {}
    for leaf in leaves:
        glyc_type = leaf.glyc_type
        d = leaf.oxygen_req_dist
        dist = pylab.array([float(d[i]) for i in order])
        dist = dist / float(sum(dist))
        matrix[glyc_type] = matrix.setdefault(glyc_type,
                                              pylab.array([0.0 for _ in order])) + dist
    
    # Normalize on the columns
    
    glyc_types = sorted(matrix.keys())
    bigram = pylab.vstack([matrix[g] for g in glyc_types])
    oxy_sums = pylab.sum(bigram, 0)
    glyc_sums = pylab.sum(bigram, 1)
    total = pylab.sum(bigram)
    int_total = int(pylab.ceil(total))
    oxy_prob = oxy_sums / float(total)
    glyc_prob = glyc_sums / float(total)
    joint_prob = pylab.outer(glyc_prob, oxy_prob)
    
    p_vals = pylab.zeros(joint_prob.shape)
    for i, glyc_type in enumerate(glyc_types):
        for j, oxygen_req in enumerate(order):
            my_pvals = []
            for k in range(int(pylab.ceil(bigram[i,j])), int_total+1):
                p = joint_prob[i,j]
                q = 1.0 - p
                choices = scipy.comb(total, k)
                pval = math.log(choices)
                pval += (k*math.log(p))
                pval += ((total-k)*math.log(q))
                my_pvals.append(math.exp(pval))
            
            p_vals[i,j] = sum(my_pvals)
            print glyc_type, oxygen_req, p_vals[i,j]
    
    
    normalized_bigram = bigram / oxy_sums
    indices = pylab.arange(len(order))
    current_bottom = pylab.zeros(len(order))
    for i, glyc in enumerate(glyc_types):
        heights = normalized_bigram[i,:]
        pylab.bar(indices, heights, color=colors[glyc],
                  bottom=current_bottom, label=glyc, width=0.5)
        current_bottom += heights
    
    # Title, ticks, labels
    #pylab.xlabel(self.projected_data.ind, fontsize='large')

    size_12 = FontProperties(size=12)
    pylab.xticks(indices + 0.25, order[:], fontproperties=size_12)
    pylab.legend()
    #axes.yaxis.set_major_locator(pylab.NullLocator())
    pylab.show()

    
    
    
if __name__ == '__main__':
    Main()