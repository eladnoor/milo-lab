#!/usr/bin/python


import pylab
import sys
import csv
from proteomics import util
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
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.genes_filename and options.counts_filename
    
    
    print 'Reading genes list from', options.genes_filename
    gene_ids = util.ReadProteinIDs(options.genes_filename)
    
    print 'Reading protein data from', options.counts_filename
    gene_counts = util.ReadProteinCounts(options.counts_filename)
    total = sum(gene_counts.values())
    print 'total count', total
    
    counts = []
    labels = []
    for id, name, count in util.ExtractCounts(gene_counts, gene_ids):
        counts.append(count)
        labels.append(name)
        
    counts = pylab.array(counts)
    pcts = counts * 100 / total
    pct = sum(pcts)
    print 'Category makes up %.2f%% of total protein' % pct
    
    fig1 = pylab.figure(0)
    pylab.hist(pcts, figure=fig1)
    
    fig2 = pylab.figure(1)
    pie_labels = ['%s %.2f%%' % t for t in zip(labels, pcts)]
    pylab.pie(counts, labels=pie_labels)
    pylab.show()
    
    
if __name__ == '__main__':
    Main()
