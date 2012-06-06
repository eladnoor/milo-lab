#!/usr/bin/python


import pylab
import sys
import csv
from proteomics import util
from toolbox.color import ColorMap
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-g", "--genes_filename",
                          dest="genes_filename",
                          help="input CSV of genes")
    opt_parser.add_option("-c", "--counts_filename",
                          dest="counts_filename",
                          help="input CSV of counts")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="input CSV of counts")
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.genes_filename and options.counts_filename
    assert options.output_filename
    
    print 'Reading genes list from', options.genes_filename
    gene_ids = util.ReadProteinIDs(options.genes_filename)
    
    print 'Reading protein data from', options.counts_filename
    gene_counts = util.ReadProteinCounts(options.counts_filename)
    total = sum(gene_counts.values())
    print 'total count', total
    
    counts = {}
    for gene_id, name, count in util.ExtractCounts(gene_counts, gene_ids):
        counts[name] = counts.get(name, 0) + count
        
    array_counts = pylab.array(counts.values())
    pcts = array_counts * 100 / total
    pct = sum(pcts)
    print 'Category makes up %.2f%% of total protein' % pct
    
    print 'Writing output CSV file to', options.output_filename
    names = sorted(set(gene_ids.values()))
    f = open(options.output_filename, 'w')
    w = csv.writer(f)
    for n in names:
        w.writerow([n, counts.get(n, 0)])
    f.close()
    
    fig1 = pylab.figure(0)
    counts_w_names = sorted([(c, n) for n,c in gene_counts.iteritems()],
                            reverse=True)
    counts_in_set = [(i, t[0]) for i, t in enumerate(counts_w_names)
                     if t[1] in gene_ids]
    
    set_idx = [t[0] for t in counts_in_set]
    set_counts = [t[1] for t in counts_in_set]
    all_counts = [t[0] for t in counts_w_names]
    max_counts = min(500, len(all_counts))
    pylab.gca().set_yscale('log')
    pylab.bar(range(max_counts), all_counts[:max_counts],
              color='g', edgecolor='g', width=2.0, figure=fig1)
    pylab.bar(set_idx, set_counts, color='r',
              figure=fig1, width=2.0)
    pylab.xlim(0, max_counts)
    
    fig2 = pylab.figure(1)
    pie_labels = sorted(counts.keys())
    colormap = ColorMap(pie_labels)
    colors = [colormap[l] for l in pie_labels]
    all_counts = [counts[k] for k in pie_labels]
    pylab.pie(all_counts, labels=pie_labels, colors=colors)
    
    for label in pie_labels:
        print label
    for count in all_counts:
        print count * 100.0 / total
        
    pylab.show()
    
    
    
    
    
if __name__ == '__main__':
    Main()
