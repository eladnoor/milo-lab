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

import genomics.stats as mystats


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


def DictSum(ds):
    if not ds:
        return {}
    out = {}
    
    for d in ds:
        for key, val in d.iteritems():
            cur_val = out.get(key, 0)
            out[key] = cur_val + val
    return out


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
    child_oxy_reqs = [c.oxygen_req for c in children]
    child_lengths = [c.edge.length for c in children]
    virtual_count = sum(c.count for c in children)
    max_length_idx = pylab.argmax(child_lengths)
    label = children[max_length_idx].taxon.label
    merged_pathways = set.union(*child_pathways)
    merged_oxygen = DictSum(child_oxy_reqs)
    
    logging.debug('Merging 2 children with edge lengths %s',
                  child_lengths)
    
    # Remove children and update the parent
    map(parent_node.remove_child, children)
    parent_node.edge.length += child_lengths[max_length_idx]
    parent_node.pathways = merged_pathways
    parent_node.count = virtual_count
    parent_node.annotate('count')
    parent_node.oxygen_req = merged_oxygen
    parent_node.annotate('oxygen_req')
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

    print 'Finding oxygen requirements'
    org_2_oxy_req = {}
    for row in db.SelectOrganismFields(['ncbi_taxon_id', 'broad_oxygen_requirement']):
        ncbi_id, oxy_req = row
        org_2_oxy_req[ncbi_id] = oxy_req
    observed_oxygen_reqs = set(org_2_oxy_req.values())
    print 'Observed oxygen requirements', observed_oxygen_reqs
    

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
        
        oxy_req = org_2_oxy_req.get(label, None)
        l.oxygen_req = {oxy_req: 1}
        l.annotate('oxygen_req')
        
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
                'upper_ed_unique': '#7b3294',
                'nonp_ed_unique': '#868800'}
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
                color = colormap.get(name, 'default')
                name_2_count[name] = name_2_count.get(name, 0) + 1
                v.append(1)
            else:
                v.append(0)
            w.writerow([taxon_id, color, value])
        path_vectors[name] = pylab.array(v)    
        f.close()
    
    ed_vector  = path_vectors['upper_ed_unique']
    emp_vector = path_vectors['upper_emp_unique']
    r, p_val = stats.pearsonr(ed_vector, emp_vector)
    print 'Pearson correlation coefficient (r)', r
    print 'R^2', r**2
    print 'p-value', p_val


    mat_shape = (len(observed_oxygen_reqs), 4)
    count_mat = np.matrix(np.zeros(mat_shape))
    prob_mat  = np.matrix(np.zeros(mat_shape))
    oxy_req_idx_map = dict((i, k) for i, k in enumerate(sorted(observed_oxygen_reqs)))
    column_labels = {0: 'None',
                     1: 'EMP Only',
                     2: 'ED Only',
                     3: 'Both'}
    
    for leaf_idx, leaf in enumerate(leaves):
        d = leaf.annotations()
        oxygen_reqs = d.get('oxygen_req', (None, None))[0]
        total_count = float(sum(oxygen_reqs.values()))
        
        ed_presence  = ed_vector[leaf_idx]
        emp_presence = emp_vector[leaf_idx]
        for i, oxy_req in oxy_req_idx_map.iteritems():
            oxy_frac = oxygen_reqs.get(oxy_req, 0) / total_count
            if ed_presence and emp_presence:
                count_mat[i, 3] += oxy_frac
            elif ed_presence:
                count_mat[i, 2] += oxy_frac
            elif emp_presence:
                count_mat[i, 1] += oxy_frac
            else:
                count_mat[i, 0] += oxy_frac
    
    total_samples = float(np.sum(count_mat))
    ed_samples    = float(np.sum(ed_vector))
    emp_samples   = float(np.sum(emp_vector))
    ed_prob       = ed_samples / total_samples
    emp_prob      = emp_samples / total_samples
    for i, oxy_req in oxy_req_idx_map.iteritems():
        oxy_req_count  = float(np.sum(count_mat[i,:]))
        oxy_req_prob   = oxy_req_count / total_samples
        # Probability of neither pathway
        prob_mat[i, 0] = (1-ed_prob) * (1-emp_prob) * oxy_req_prob
        # Probability of EMP only
        prob_mat[i, 1] = (1-ed_prob) * (emp_prob) * oxy_req_prob
        # Probability of ED only
        prob_mat[i, 2] = (ed_prob) * (1-emp_prob) * oxy_req_prob
        # Probability of both
        prob_mat[i, 3] = (ed_prob) * (emp_prob) * oxy_req_prob
    
    p_vals = mystats.CalcPValues(count_mat, prob_mat)
    print 'Counts'
    print count_mat
    print 'Probabilities assuming random sampling'
    print prob_mat
    print 'Mappings'
    print oxy_req_idx_map
    print column_labels
    print 'P-values'
    print p_vals
    
    # Plot the p-values
    pylab.figure()
    xs = sorted(column_labels.keys())
    xticks = [column_labels[i] for i in xs]
    ys = sorted(oxy_req_idx_map.keys())
    yticks = [oxy_req_idx_map[j] for j in ys]
    
    sigs = p_vals < 0.05
    super_sigs = p_vals < 0.001
    rows, cols = sigs.shape
    for i in xrange(rows):
        for j in xrange(cols):
            if super_sigs[i,j]:
                print oxy_req_idx_map[i], 'x', column_labels[j],
                print '**', p_vals[i,j]
                pylab.text(j, i, '**', color='w')
            elif sigs[i,j]:
                print oxy_req_idx_map[i], 'x', column_labels[j],
                print '*', p_vals[i,j]
                pylab.text(j, i, '*', color='w')
    
    
    pylab.imshow(p_vals, interpolation='nearest')
    pylab.xticks(xs, xticks)
    pylab.yticks(ys, yticks)
    pylab.colorbar()
    
    
    # Plot the bar plot breakdown.
    
    restricted_counts = np.matrix(np.zeros((3,3)))
    allowed_genotypes = ['ED Only', 'Both', 'EMP Only']
    allowed_phenotypes = ['anaerobe', 'facultative', 'aerobe']
    
    for j, genotype in column_labels.iteritems():
        if genotype not in allowed_genotypes:
            continue
        
        new_j = allowed_genotypes.index(genotype)
        
        for i, phenotype in oxy_req_idx_map.iteritems():
            if phenotype not in allowed_phenotypes:
                continue
            
            new_i = allowed_phenotypes.index(phenotype)
            restricted_counts[new_i, new_j] = count_mat[i,j]
    
    pcts_matrix = restricted_counts / np.sum(restricted_counts, 1)
    print 'Phenotypes (rows)'
    print allowed_phenotypes
    print 'Genotypes (cols)'
    print allowed_genotypes
    print 'Counts of interesting categories'
    print restricted_counts
    print 'PCTs including interesting categories'
    print pcts_matrix * 100.0
    print 'Effective number of organisms', np.sum(restricted_counts)
    
    colors = ['#37DD6F', '#FF5D40', '#4186D3']
    pylab.figure()
    current_bottom = pylab.zeros(3)
    rows, cols = pcts_matrix.shape
    for j in xrange(cols):
        heights = np.array(pcts_matrix[:,j].flat)
        xs = np.arange(heights.size)
        pylab.bar(xs, heights, bottom=current_bottom,
                  color=colors[j], edgecolor='w',
                  label=allowed_genotypes[j],
                  align='center')
        current_bottom += heights
    xs = pylab.arange(3) + 0.5
    pylab.xticks(xs, allowed_phenotypes)
    pylab.legend()
    pylab.show()
    
    
    colormap = {'Organic': '#ff0000',
                'Inorganic': '#00ff00',
                'aerobe': '#2861e4',
                'anaerobe': '#e44228',
                'facultative': '#e4c328'}
    fname = 'trophism.csv'
    f = open(fname, 'w')
    w = csv.writer(f)
    for l in leaves:
        label = l.taxon.label
        cat = db.NCBI2EnergyCategory(label)
        color = colormap.get(cat, default_color)
        w.writerow([label, color, cat])
    f.close()

    fname = 'oxy_req.csv'
    f = open(fname, 'w')
    w = csv.writer(f)
    for l in leaves:
        label = l.taxon.label
        cat = db.NCBI2BroadOxygenReq(label)
        color = colormap.get(cat, default_color)
        w.writerow([label, color, cat])
    f.close()
    
    nleaves = len(leaves)
    for name, count in name_2_count.iteritems():
        pct = 100.0 * float(count) / float(nleaves)
        print '%.2f%% (%d of %d) have pathway %s' % (pct, count, nleaves,
                                                    str(name))
    
    v_or_accumulator  = pylab.zeros(nleaves)
    v_and_accumulator = pylab.ones(nleaves)
    for v in path_vectors.values():
        v_or_accumulator  = np.logical_or(v_or_accumulator, v)
        v_and_accumulator = np.logical_and(v_and_accumulator, v)
    total_w_any = pylab.where(v_or_accumulator == True)[0].size
    total_w_all = pylab.where(v_and_accumulator == True)[0].size
    any_pct = 100.0 * float(total_w_any) / float(nleaves)
    all_pct = 100.0 * float(total_w_all) / float(nleaves)
    print '%.2f%% (%d of %d) have some pathway' % (any_pct,
                                                   total_w_any,
                                                   nleaves)
    print '%.2f%% (%d of %d) have all pathways' % (all_pct,
                                                   total_w_all,
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