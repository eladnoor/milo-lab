#!/usr/bin/python

import pylab
import sys

from genomics.monte_carlo_tester import MonteCarloTester
from genomics.raw_data import RawData
from genomics.row_filterer import RowFilterer
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          help="input CSV")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="output CSV")
    opt_parser.add_option("-a", "--first_col",
                          dest="first_col",
                          help="First column in 2-way histogram")
    opt_parser.add_option("-b", "--second_col",
                          dest="second_col",
                          help="Second column in 2-way histogram")
    opt_parser.add_option("-c", "--third_col",
                          dest="third_col",
                          help="Third column in 2-way histogram")
    opt_parser.add_option("-v", "--second_col_key",
                          dest="second_col_key",
                          help="First column in 2-way histogram")
    opt_parser.add_option("-w", "--third_col_key",
                          dest="third_col_key10000",
                          help="First column in 2-way histogram")    
    opt_parser.add_option("-f", "--filter_values",
                          dest="filter_values",
                          help="Value to ignore.")
    opt_parser.add_option("-t", "--title",
                          dest="title",
                          help="Figure Title")
    opt_parser.add_option("-x", "--filter_cols",
                          dest="filter_cols",
                          help="Columns to filter on, comma separated")
    opt_parser.add_option("-y", "--filter_cols_vals",
                          dest="filter_cols_vals",
                          help="Values to keep, one per column")
    opt_parser.add_option("-n", "--num_trials",
                          dest="num_trials",
                          default=1000,
                          type="int",
                          help="Number of monte-carlo samples.")
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename
    assert options.first_col and options.second_col
    print 'Reading species list from', options.input_filename

    # Read options
    first_col, second_col = options.first_col, options.second_col
    third_col = options.third_col
    filter_values = []
    if options.filter_values:
        filter_values = set(map(str.strip, options.filter_values.split(',')))

    # Set up row filter    
    row_filterer = None
    if options.filter_cols:
        assert options.filter_cols_vals
        filter_cols = map(str.strip, options.filter_cols.split(','))
        filter_cols_vals = map(str.strip, options.filter_cols_vals.split(','))
        row_filterer = RowFilterer(filter_cols, filter_cols_vals)

    # Set up dependent columns to read.
    dependent_cols = [second_col]
    if third_col:
        dependent_cols.append(third_col)
    
    # Read data and compute histogram.
    data = RawData.FromCsvFile(options.input_filename,
                               first_col,
                               dependent_cols,
                               row_filterer=row_filterer)
    projected = data.Project()
    histogram = projected.MakeHistogram(filter_values=filter_values)
    
    print len(histogram.pairs), 'filtered examples'
    
    # Unfiltered overall distribution.
    dep_counts = projected.DepCounts()
    
    deps = sorted(dep_counts.keys())
    dcounts = [dep_counts[l] for l in deps]
    dlabels = [l or 'None given' for l in deps]
    fig4 = pylab.figure(4)
    pylab.pie(dcounts, labels=dlabels, autopct='%.1f%%')
    
    # Filtered overall distribution.
    ind_counts = histogram.FilteredIndCounts()
    dep_counts = histogram.FilteredDepCounts()
    
    inds = sorted(ind_counts.keys())
    icounts = [ind_counts[l] for l in inds]
    ilabels = [l or 'None given' for l in inds]
    fig0 = pylab.figure(0)
    pylab.pie(icounts, labels=ilabels, autopct='%.1f%%')

    deps = sorted(dep_counts.keys())
    dcounts = [dep_counts[l] for l in deps]
    dlabels = [l or 'None given' for l in deps]
    fig1 = pylab.figure(1)
    pylab.pie(dcounts, labels=dlabels, autopct='%.1f%%')
    
    # Calculate and print p-values
    sampler = MonteCarloTester(histogram)
    sampler.Test(options.num_trials)
    fig2 = pylab.figure(2)
    sampler.HeatMap(fig2)
    
    # Bar plot histogram
    fig3 = pylab.figure(3)
    axes = pylab.axes()
    fig3.add_axes(axes)
    title = options.title or '%s vs. %s' % (options.first_col, options.second_col)
    pylab.title(title, figure=fig3)
    
    histogram.BarPlot(axes)
    pylab.legend()
    pylab.show()
    



if __name__ == '__main__':
    Main()
