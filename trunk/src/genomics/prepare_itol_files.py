#!/usr/bin/python

import csv
import dendropy
import logging
import pylab
import sys

import numpy as np

from collections import Counter
from genomics import genome_db
from genomics import pathway
from toolbox import color
from toolbox import util
from optparse import OptionParser
from scipy import stats


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-t", "--tree_filename",
                          dest="tree_filename",
                          help="Phylogenetic tree filename")
    opt_parser.add_option("-f", "--tree_format",
                          dest="tree_format",
                          default="newick",
                          help="Phylogenetic tree format")
    opt_parser.add_option("-p", "--pathways_filename",
                          dest="pathways_filename",
                          default='../data/genomics/glycolysis_pathways_unique.json',
                          help="Input pathways JSON file")
    opt_parser.add_option("-d", "--genome_db_filename",
                          dest="genome_db_filename",
                          default='../res/genomes.sqlite',
                          help="Genome database filename")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="output annotated and pruned tree")
    opt_parser.add_option("-e", "--edge_length_threshold",
                          type="float",
                          default=0.0,
                          dest="edge_length_threshold",
                          help="Edge length threshold.")
    opt_parser.add_option("-r", "--only_heterotrophs",
                          dest="only_heterotrophs",
                          action="store_true",
                          default=False,
                          help="Whether to filter out non-heterotrophs.")
    return opt_parser



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
    child_pathways = [c.pathways for c in children]
    child_lengths = [c.edge.length for c in children]
    virtual_count = sum(c.count for c in children)
    max_length_idx = pylab.argmax(child_lengths)
    label = children[max_length_idx].taxon.label
    merged_pathways = set.union(*child_pathways)
    
    logging.debug('Merging 2 children with edge lengths %s',
                  child_lengths)
    
    # Remove children and update the parent
    map(parent_node.remove_child, children)
    parent_node.edge.length += child_lengths[max_length_idx]
    parent_node.pathways = merged_pathways
    parent_node.count = virtual_count
    parent_node.annotate('count')
    for pname in parent_node.pathways:
        setattr(parent_node, pname, True)
        parent_node.annotate(pname)
    
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
    assert options.pathways_filename and options.genome_db_filename
    assert options.output_filename
    
    print 'Reading tree from', options.tree_filename
    tree = dendropy.Tree()
    tree.read_from_path(options.tree_filename, options.tree_format,
                        extract_comment_metadata=True)
    
    leaves = tree.leaf_nodes()
    print 'Tree has %d leaf nodes' % len(leaves)
    
    pathways = pathway.LoadPathways(options.pathways_filename)
    db = genome_db.GenomeDB(options.genome_db_filename)
    
    pathway_names = [util.slugify(p.name) for p in pathways]
    org_2_pathways = {}
    for path in pathways:
        org_counts = Counter()
        for enz_set in path.enzyme_sets:
            orgs_w_enz = set()
            for ec in enz_set:
                orgs_w_enz.update(list(db.OrganismsForEC(ec)))                
            org_counts += Counter(orgs_w_enz)
            
        orgs_w_pathway = [o for o,c in org_counts.iteritems()
                          if c == len(path.enzyme_sets)]
        orgs_w_pathway = filter(None, map(db.KEGG2NCBI, orgs_w_pathway))
        
        for org in orgs_w_pathway:
            org_2_pathways.setdefault(org, set()).add(util.slugify(path.name))

    # Find the organisms that have pathway tags.         
    all_labels = set([l.taxon.label for l in leaves])
    pathway_orgs = set(org_2_pathways.keys())
    intersect = all_labels.intersection(pathway_orgs)
    print len(intersect), 'pathway orgs found'
    print len(pathway_orgs) - len(intersect), 'pathway orgs not found'
    
    # Find organisms that are heterotrophs
    if options.only_heterotrophs:
        print 'Pruning non-heterotrophs'
        q = db.db.Execute('SELECT ncbi_taxon_id, energy_category from organisms')
        ncbi_to_keep = set()
        for row in q:
            ncbi_id, energy_cat = row
            if energy_cat and energy_cat.lower() == 'organic':
                ncbi_to_keep.add(ncbi_id.strip())
        tree.retain_taxa_with_labels(ncbi_to_keep)
    
        leaves = tree.leaf_nodes()
        print 'Tree now contains', len(leaves), 'leaves'
    
    lengths = []
    for e in tree.leaf_edge_iter():
        lengths.append(e.length)
    lengths = pylab.array(lengths)
    below_thresh = pylab.find(lengths < options.edge_length_threshold).size
    pct_below = 100.0 * float(below_thresh) / float(len(lengths))
    
    print 'Median length', pylab.median(lengths)
    print 'Mean length', pylab.mean(lengths)
    print pct_below, '% below threshold'
    
    print 'Pruning leaves'            
    for l in tree.leaf_nodes():
        label = l.taxon.label
        pathways = org_2_pathways.get(label, set())
        l.pathways = pathways
        for pname in pathways:
            setattr(l, pname, True)
            l.annotate(pname)
        l.count = 1
        l.annotate('count')

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

    colormap = {'upper_emp_unique': '#008837',
                'upper_ed_unique': '#7b3294'}
    default_color = '#c0c3c7'
    
    name_2_count = {}
    path_vectors = {}
    for name in pathway_names:
        fname = util.slugify(name) + '.csv'
        f = open(fname, 'w')
        w = csv.writer(f)
        v = []
        for leaf in leaves:
            taxon_id = leaf.taxon.label
            pathways = leaf.pathways
            d = leaf.annotations()
            value = d.get(name, (False, None))[0]
            color = default_color
            if value is True:
                color = colormap[name]
                name_2_count[name] = name_2_count.get(name, 0) + 1
                v.append(1)
            else:
                v.append(0)
            w.writerow([taxon_id, color, value])
        path_vectors[name] = pylab.array(v)    
        f.close()
    
    r, p_val = stats.pearsonr(path_vectors['upper_emp_unique'],
                              path_vectors['upper_ed_unique'])
    print 'Pearson correlation coefficient (r)', r
    print 'R^2', r**2
    print 'p-value', p_val

    colormap = {'Organic': '#ff0000',
                'Inorganic': '#00ff00',
                None: '#0000ff'}
    fname = 'trophism.csv'
    f = open(fname, 'w')
    w = csv.writer(f)
    for l in leaves:
        label = l.taxon.label
        cat = db.NCBI2EnergyCategory(label)
        color = colormap.get(cat, '#0000ff')
        w.writerow([label, color, cat])
    f.close()
    
    nleaves = len(leaves)
    for name, count in name_2_count.iteritems():
        pct = 100.0 * float(count) / float(nleaves)
        print '%.2f%% (%d of %d) have pathway %s' % (pct, count, nleaves,
                                                    str(name))
    
    v_accumulator = pylab.zeros(nleaves)
    for v in path_vectors.values():
        v_accumulator = np.logical_or(v_accumulator, v)
    total_w_any = pylab.where(v_accumulator == True)[0].size
    any_pct = 100.0 * float(total_w_any) / float(nleaves)
    print '%.2f%% (%d of %d) have some pathway' % (any_pct,
                                                   total_w_any,
                                                   nleaves)
    
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