#!/usr/bin/python


import pylab
import sys
import cvxmod

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
    opt_parser.add_option("-l", "--a_label",
                          dest="a_label",
                          default="Protein levels A",
                          help="label for levels A")
    opt_parser.add_option("-b", "--protein_levels_b",
                          dest="protein_levels_b",
                          help="input CSV of second protein levels")
    opt_parser.add_option("-m", "--b_label",
                          dest="b_label",
                          default="Protein levels B",
                          help="label for levels B")
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.genes_filename
    assert options.protein_levels_a and options.protein_levels_b
    
    
    print 'Reading genes list from', options.genes_filename
    gene_ids = util.ReadProteinIDs(options.genes_filename)
    
    print 'Reading protein data A from', options.protein_levels_a
    gene_counts_a = util.ReadProteinCounts(options.protein_levels_a)
    print 'Reading protein data B from', options.protein_levels_b
    gene_counts_b = util.ReadProteinCounts(options.protein_levels_b)

    my_counts_a = dict((id, (count, name)) for id, name, count in
                       util.ExtractCounts(gene_counts_a, gene_ids))
    my_counts_b = dict((id, (count, name)) for id, name, count in
                       util.ExtractCounts(gene_counts_b, gene_ids))
        
    overlap_ids = set(my_counts_a.keys()).intersection(my_counts_b.keys())
    x = pylab.matrix([my_counts_a[id][0] for id in overlap_ids])
    y = pylab.matrix([my_counts_b[id][0] for id in overlap_ids])
    labels = [my_counts_b[id][1] for id in overlap_ids]
    
    xlog = pylab.log10(x)
    ylog = pylab.log10(y)
    a = cvxmod.optvar('a', 1)
    mx = cvxmod.matrix(xlog.T)
    my = cvxmod.matrix(ylog.T)
    
    p = cvxmod.problem(cvxmod.minimize(cvxmod.norm2(my - a - mx)))
    p.solve(quiet=True)
    offset = cvxmod.value(a)
    lin_factor = 10**offset
    lin_label = 'Y = %.2g*X' % lin_factor
    log_label = 'log10(Y) = %.2g + log10(X)' % offset
    
    f1 = pylab.figure(0)
    pylab.title('Linear scale')
    xylim = max([x.max(), y.max()]) + 5000
    linxs = pylab.arange(0.0, xylim, 0.1)
    linys = linxs * lin_factor
    pylab.plot(x.tolist()[0], y.tolist()[0], 'g.', label='Protein Data')
    pylab.plot(linxs, linys, 'b-', label=lin_label)
    for x_val, y_val, label in zip(x.tolist()[0], y.tolist()[0], labels):
        pylab.text(x_val, y_val, label, fontsize=8)

    pylab.xlabel(options.a_label)
    pylab.ylabel(options.b_label)
    pylab.legend()
    pylab.xlim((0.0, xylim))
    pylab.ylim((0.0, xylim))
    
    f2 = pylab.figure(1)
    pylab.title('Log10 scale')
    xylim = max([xlog.max(), ylog.max()]) + 1.0    
    pylab.plot(xlog.tolist()[0], ylog.tolist()[0], 'g.', label='Log10 Protein Data')
    linxs = pylab.arange(0.0, xylim, 0.1)
    linys = linxs + offset
    pylab.plot(linxs, linys, 'b-', label=log_label)
    
    for x_val, y_val, label in zip(xlog.tolist()[0], ylog.tolist()[0], labels):
        pylab.text(x_val, y_val, label, fontsize=8)
    
    pylab.xlabel(options.a_label + ' (log10)')
    pylab.ylabel(options.b_label + ' (log10)')
    pylab.legend()
    pylab.xlim((0.0, xylim))
    pylab.ylim((0.0, xylim))
    
    
    
    pylab.show()
        
    
    
if __name__ == '__main__':
    Main()
