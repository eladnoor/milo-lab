#!/usr/bin/python

import itertools
import pylab

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
    
    def GetDepKeys(self):
        dep_keys = set()
        for count_d in self.counts.values():
            dep_keys.update(count_d.keys())
        return dep_keys
    
    def BarPlot(self, axes):
        weight_array, count_array = [], []
        ind_keys = self.ind_keys
        dep_keys = sorted(self.dep_keys)
        for ind in ind_keys:
            weight_array.append(self.weights.get(ind, {}))
            count_array.append(self.counts.get(ind, {}))        

        colormap = ColorMap(dep_keys)
        indices = pylab.arange(len(ind_keys))
        current_bottom = pylab.zeros(len(ind_keys))
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
        labels = [a or "None given" for a in ind_keys]
        pylab.xticks(indices + 0.25, labels, fontproperties=size_12)
        axes.yaxis.set_major_locator(pylab.NullLocator())