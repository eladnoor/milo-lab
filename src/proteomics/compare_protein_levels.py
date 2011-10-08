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
    opt_parser.add_option("-a", "--protein_levels_a",
                          dest="protein_levels_a",
                          help="input CSV of protein levels")
    opt_parser.add_option("-b", "--protein_levels_b",
                          dest="protein_levels_b",
                          help="input CSV of second protein levels")
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.genes_filename
    assert options.protein_levels_a and options.protein_levels_b
    
    
    print 'Reading genes list from', options.genes_filename
    gene_ids = util.ReadProteinIDs(options.genes_filename)
    
    print 'Reading protein data from', options.protein_levels_a
    gene_counts_a = util.ReadProteinCounts(options.protein_levels_a)
    total_a = sum(gene_counts_a.values())
    print 'Reading protein data from', options.protein_levels_b
    gene_counts_b = util.ReadProteinCounts(options.protein_levels_b)
    total_b = sum(gene_counts_b.values())

    my_counts_a = dict((id, (count, name)) for id, name, count in
                       util.ExtractCounts(gene_counts_a, gene_ids))
    my_counts_b = dict((id, (count, name)) for id, name, count in
                       util.ExtractCounts(gene_counts_b, gene_ids))
    
    overlap_ids = set(my_counts_a.keys()).intersection(my_counts_b.keys())
    x = pylab.matrix([my_counts_a[id][0] for id in overlap_ids])
    y = pylab.matrix([my_counts_b[id][0] for id in overlap_ids])
    labels = [my_counts_b[id][1] for id in overlap_ids]
    x = pylab.log(x)
    y = pylab.log(y)
    
    print x
    print y
    mx = pylab.vstack([x, pylab.ones(x.shape)]).T
    my = pylab.matrix(y).T
    print mx.shape
    print my.shape
    a = pylab.lstsq(mx, my)[0]
    slope = a[0, 0]
    offset = a[1, 0]
    print slope, offset

    xylim = max([x.max(), y.max()])
    pylab.xlim((0.0, xylim))
    pylab.ylim((0.0, xylim))
    
    pylab.plot(x, y, 'g.')
    linxs = pylab.arange(0.0, xylim, 0.1)
    linys = slope * linxs + offset
    print linxs, linys
    pylab.plot(linxs, linys, 'b-')
    for x_val, y_val, label in zip(x, y, labels):
        pylab.text(x_val+0.1, y_val+0.1, label, fontsize=8)
    

    pylab.show()
        
    
    
if __name__ == '__main__':
    Main()
