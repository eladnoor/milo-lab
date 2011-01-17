import pylab
from scipy.stats.stats import ranksums
import matplotlib

def cdf(v, label, style=None):
    new_v = []
    for x in v:
        if (x != None):
            new_v.append(x)
    l = len(new_v)
    if (l < 2):
        pylab.plot([0, 0], [0, 0], label=label)
    else:
        new_v.sort()
        
        if (style == None):
            pylab.plot(new_v, [float(x)/(float(l)-1) for x in range(l)], label=label)
        else:
            pylab.plot(new_v, [float(x)/(float(l)-1) for x in range(l)], style, label=label)

def query_ranksum(cursor, query1, query2):
    x_list1 = []
    x_list2 = []
    
    cursor.execute(query1)
    for row in cursor:
        x_list1.append(row[0])
    
    cursor.execute(query2)
    for row in cursor:
        x_list2.append(row[0])
    
    (z_statistic, p_value) = ranksums(x_list1, x_list2)
    return p_value

def plot_cdf(cursor, query, prefix='cdf', xlog=False, xlimits=None, xlabel='', fontsize=10, title=''):
    """
        Plots the CDFs of a series of values according to the first column of the query.
        The second column is a filter label over which to draw the distribution 
    """
    styles = ['b-', 'r-', 'g-', 'c-', 'y-', 'm-', 'b--', 'r--', 'g--', 'c--', 'y--', 'm--']
    cursor.execute(query)
    x_lists = {}
    total_list = []
    for row in cursor:
        (x, l) = row
        l = str(l)
        if (not l in x_lists):
            x_lists[l] = []
        x_lists[l].append(x)
        total_list.append(x)
    
    pylab.figure()
    pylab.hold(True)
    cdf(total_list, 'total', style='k:')
    counter = 0
    for l in sorted(x_lists.keys()):
        cdf(x_lists[l], l + " (%d)" % len(x_lists[l]), style=styles[counter])
        counter += 1

    if (xlog):
        pylab.xscale('log')
    
    if (xlimits != None):
        pylab.xlim(xlimits)
    
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel('CDF (%)')
    legendfont = matplotlib.font_manager.FontProperties(size=fontsize)
    pylab.legend(loc='lower right', prop=legendfont)
    pylab.savefig('../res/%s.pdf' % prefix, format='pdf')
    pylab.hold(False)
            
def plot_xy(cursor, query, prefix='plotxy', color='b', marker='.', xlog=False, ylog=False, xlabel='', ylabel='', title=''):
    """
        Executes the 'query' which should return two numerical columns.
    """
    cursor.execute(query)
    x_list = []
    y_list = []
    for row in cursor:
        (x, y) = row
        x_list.append(x)
        y_list.append(y)
        
    pylab.figure()
    pylab.hold(True)
    pylab.plot(x_list, y_list, color=color, marker=marker, linestyle='None', markersize=5)
    if (xlog):
        pylab.xscale('log')
    if (ylog):
        pylab.yscale('log')
        
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    pylab.savefig('../res/%s.pdf' % prefix, format='pdf')
    pylab.hold(False)

def plot_scatter(cursor, query, prefix='scatter', marker='.', xlog=False, ylog=False, xlabel='', ylabel='', fontsize=10, title=''):
    """
        Executes the 'query' which should return two numerical columns and an additional label column
    """
    cursor.execute(query)
    x_lists = {}
    y_lists = {}
    styles = ['rx','gx','bx','yx','mx','kx','cx','r.','g.','b.','y.','m.','k.','c.']
    for row in cursor:
        (x, y, l) = row
        if (not l in x_lists):
            x_lists[l] = []
            y_lists[l] = []
        x_lists[l].append(x)
        y_lists[l].append(y)
    
    pylab.figure()
    pylab.hold(True)
    counter = 0
    for l in sorted(x_lists.keys()):
        pylab.plot(x_lists[l], y_lists[l], styles[counter % len(styles)], label=str(l), markersize=5)
        counter += 1

    if (xlog):
        pylab.xscale('log')
    if (ylog):
        pylab.yscale('log')
        
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    legendfont = matplotlib.font_manager.FontProperties(size=fontsize)
    pylab.legend(loc='lower right', prop=legendfont)
    pylab.savefig('../res/%s.pdf' % prefix, format='pdf')
    pylab.hold(False)
    
def plot_bar(cursor, query, prefix='plotbar', xlabel='', ylabel=''):
    """
        Executes the 'query' which should return a textual column followed by a numerical column.
    """
    cursor.execute(query)
    x_list = []
    y_list = []
    for row in cursor:
        (x, y) = row
        x_list.append(x)
        y_list.append(y)
        
    pylab.figure()
    pylab.hold(True)
    ind = pylab.arange(len(x_list))
    width = 0.35
    pylab.bar(ind, y_list)
    pylab.xticks(ind+width/2., tuple(x_list))
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    pylab.savefig('../res/%s.pdf' % prefix, format='pdf')
    #pylab.show()
    pylab.hold(False)

def median(v):
    l = len(v)
    if (l == 0):
        return None
    v.sort()
    if (l % 2 == 1):
        return v[(l-1)/2]
    else:
        return (v[l/2 - 1] + v[l/2]) * 0.5

def median_group_by(cursor, query, group_by_cols, median_cols):
    cursor.execute(query)
    query_map = {}
    for row in cursor:
        key = tuple([row[c] for c in group_by_cols])
        value = tuple([row[c] for c in median_cols])
        if (not key in query_map):
            query_map[key] = [value]
        else:
            query_map[key].append(value)
    
    result_table = []
    for (key, values) in sorted(query_map.iteritems()):
        res = list(key)
        for i in range(len(median_cols)):
            v = []
            for j in range(len(values)):
                if (values[j][i] != None):
                    v.append(values[j][i])
            res.append(median(v))
        result_table.append(res)
    
    return result_table

def count_unique(cursor, query):
    counters = {}
    counters['total'] = 0
    for row in cursor.execute(query):
        counters[tuple(row)] = counters.get(tuple(row), 0) + 1
        counters['total'] += 1
    return counters
