import pylab
import sys
import csv
from toolbox.color import ColorMap
from optparse import OptionParser

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          help="input CSV of counts")
    return opt_parser

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename
    
    f = open(options.input_filename)
    r = csv.DictReader(f)
    orgs = r.fieldnames[1:]
    names = []
    data = []
    for i, row in enumerate(r):
        names.append(row.get('Enzyme Function'))
        a = pylab.array([float(row[k]) for k in orgs])
        data.append(a)

    

    m = pylab.array(data).T
    rows, cols = m.shape
    cur_bottom = pylab.zeros(rows)
    left = range(rows)
    colors = ColorMap(range(cols))
    for i in xrange(cols):
        heights = m[:,i]
        pylab.bar(left, heights, bottom=cur_bottom, color=colors[i],
                  width=0.5, align='center', label=names[i])
        cur_bottom += heights
    pylab.ylim((0.0, cur_bottom.max()))
    pcts = pylab.sum(m, 1)
    ticks = ['%s \n%.1f%%' % (o,p) for o,p in zip(orgs,pcts)]
    pylab.xticks(pylab.arange(rows), ticks, ha='center')
    pylab.legend()
    pylab.show()
    
    
if __name__ == '__main__':
    Main()