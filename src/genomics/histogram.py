#!/usr/bin/python

import itertools
import scipy
import pylab
import math

from genomics.colormap import ColorMap
from matplotlib.font_manager import FontProperties


class Histogram(object):
    
    def __init__(self, projected_data, filter_values=None):
        self.projected_data = projected_data
        self.raw_data = projected_data.raw_data
        self.filter_values = filter_values
        self.pairs = list(self.projected_data.Iterate(filter_values=self.filter_values))
        self.counts = self.MakeCounts()
        self.weights = self.MakeWeights()
        self.dep_keys = self.GetDepKeys()
        self.ind_keys = sorted(self.counts.keys())
    
    def AllPossiblePairs(self):
        return itertools.product(self.ind_keys, self.dep_keys)
    
    def FilteredIndCounts(self):
        """Post-filtered counts."""
        c = {}
        for ind, unused_dep in self.pairs:
            c[ind] = c.get(ind, 0) + 1
        return c

    def FilteredDepCounts(self):
        """Post-filtered counts."""
        c = {}
        for unused_ind, dep in self.pairs:
            c[dep] = c.get(dep, 0) + 1
        return c
    
    def GetWeight(self, ind, dep):
        return self.weights.get(ind, {}).get(dep, 0.0)
    
    def MakeCounts(self):
        counts = {}
        for ind, dep in self.pairs:
            val_dict = counts.setdefault(ind, {})
            val_dict[dep] = val_dict.get(dep, 0) + 1
        return counts
    
    def MakeWeights(self):
        weights = {}
        for ind, count_d in self.counts.iteritems():
            total = float(sum(count_d.values()))
            weight_d = dict((k, float(v)/total) for k,v in count_d.iteritems())
            weights[ind] = weight_d
            
        return weights
    
    def CalcPValues(self):
        filtered_ind = filter(None, self.ind_keys)
        filtered_dep = filter(None, self.dep_keys)
        filtered_pairs = filter(lambda x: x[0] != None and x[1] != None,
                                self.pairs)
        m = pylab.zeros((len(filtered_ind), len(filtered_dep)))
        total = len(filtered_pairs)
        total2 = float(total**2)
        for i, ind_key in enumerate(filtered_ind):
            for j, dep_key in enumerate(filtered_dep):
                ind_counts = self.counts.get(ind_key, {})
                
                total_ind_count = float(sum(ind_counts.itervalues()))
                print ind_key, ':', total_ind_count / float(total)
                total_dep_count = float(len(filter(lambda x: x[1] == dep_key,
                                                   self.pairs)))
                print dep_key, ':', total_dep_count / float(total)
                
                prob = total_ind_count * total_dep_count / total2
                observed_successes = ind_counts.get(dep_key, 0)
                pvals = []
                for successes in xrange(observed_successes, total + 1):
                    
                    pval = math.log(scipy.comb(total, successes, exact=True))
                    pval += (successes*math.log(prob))
                    pval += ((total-successes)*math.log(1-prob))
                    pvals.append(math.exp(pval))
                pval = sum(pvals)
                print ind_key, dep_key, pval
                m[i][j] = pval
        
        return m
    
    def PValHeatMap(self, figure):
        pmat = self.CalcPValues()
        print pmat
        
        max_x, max_y = pmat.shape
        print max_x, max_y
        
        dep_labels = filter(None, self.dep_keys)
        ind_labels = filter(None, self.ind_keys)
                
        pylab.imshow(pmat, interpolation="nearest", figure=figure)
        nx = float(len(dep_labels))
        ny = float(len(ind_labels))
        pylab.xticks(pylab.arange(nx), dep_labels, rotation=20, ha='right', figure=figure)
        pylab.yticks(pylab.arange(ny), ind_labels, figure=figure)
        
        idx = zip(*(pmat < 1e-2).nonzero())
        print idx
        for y,x in idx:
            s = '*'
            if pmat[y][x] < 1e-5:
                s = '**'
            pylab.text(x, y, s,
                       color='w', weight='bold', ha='center')
        
        pylab.title('p-values', figure=figure)
        pylab.colorbar()
        pylab.axis('scaled')
        
    def GetDepKeys(self):
        dep_keys = set()
        for count_d in self.counts.values():
            dep_keys.update(count_d.keys())
        return dep_keys
    
    def BarPlot(self, axes, show_ind_none=False):
        weight_array, count_array = [], []
        labels = []
        ind_keys = self.ind_keys
        dep_keys = sorted(self.dep_keys)
        for ind in ind_keys:
            if not ind and not show_ind_none:
                continue
                
            weight_array.append(self.weights.get(ind, {}))
            count_array.append(self.counts.get(ind, {}))
            labels.append(ind or "None given")        

        colormap = ColorMap(dep_keys)
        indices = pylab.arange(len(weight_array))
        current_bottom = pylab.zeros(len(weight_array))
        for dep in dep_keys:
            heights = pylab.array([w.get(dep, 0.0) for w in weight_array])
            pylab.bar(indices, heights, color=colormap[dep],
                      bottom=current_bottom, label=dep, width=0.5)
            
            counts = [c.get(dep, 0) for c in count_array]
            for x, y, count, weight in zip(indices, current_bottom, counts, heights):
                if count == 0:
                    continue
                
                txt = '%d; %.1f%%' % (count, 100*weight)
                pylab.text(x + 0.25, y + 0.02, txt, ha='center')
                
            current_bottom += heights
        
        # Title, ticks, labels
        pylab.xlabel(self.projected_data.ind, fontsize='large')

        size_12 = FontProperties(size=12)
        pylab.xticks(indices + 0.25, labels, fontproperties=size_12)
        axes.yaxis.set_major_locator(pylab.NullLocator())