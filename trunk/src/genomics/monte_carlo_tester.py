#!/usr/bin/python

import itertools
import pylab


class MonteCarloTester(object):
    """Calculates p-values for a given histogram using monte-carlo testing."""
    def __init__(self, hist, filter_values=None):
        """Initialize.
        
        Args:
            hist: the histogram object to compute p-vals for.
        """
        self.hist = hist
        self.raw_data = hist.raw_data
        self.filter_values = hist.filter_values
        self.all_pairs = sorted(self.hist.AllPossiblePairs())
        self.all_ind = list(set([p[0] for p in self.all_pairs]))
        self.all_dep = list(set([p[1] for p in self.all_pairs]))
        self.greater_p = None
        self.lower_p = None
        
    def Test(self, n=1000):
        """Test the histogram.
        
        Args:
            n: the number of trials. Defaults to 1000.
            
        Returns:
            A 2-tuple (greater p-vals, lower p-vals) where each 
            member is a dictionary mapping a pair (independent, dependent) to
            a p-value. The first dict contains p-values representing the 
            probability of observing a greater value and the second lesser.
        """
        ge_counts, le_counts = {}, {}
        all_pairs = sorted(self.hist.AllPossiblePairs())
        sampled_weights = dict((p, []) for p in all_pairs)
        for _ in xrange(n):
            shuffled = self.raw_data.Shuffle()
            shuffled_hist = shuffled.Project().MakeHistogram(
                filter_values=self.filter_values)
            
            for pair in all_pairs:
                ind, dep = pair
                observed_weight = self.hist.GetWeight(ind, dep)
                shuffled_weight = shuffled_hist.GetWeight(ind, dep)
                sampled_weights[pair].append(shuffled_weight)
                
                if shuffled_weight >= observed_weight:
                    ge_counts[pair] = ge_counts.get(pair, 0) + 1
                if shuffled_weight <= observed_weight:
                    le_counts[pair] = le_counts.get(pair, 0) + 1
        
        total = float(n)
        self.greater_p = dict((p, float(ge_counts.get(p, 0.0))/total) for p in all_pairs)
        self.lower_p = dict((p, float(le_counts.get(p, 0.0))/total) for p in all_pairs)
        
        return self.greater_p, self.lower_p
    
    def MakeMatrix(self, p_dict):
        m = pylab.zeros((len(self.all_ind), len(self.all_dep)))
        for i, ind in enumerate(self.all_ind):
            for j, dep in enumerate(self.all_dep):
                m[i][j] = p_dict.get((ind, dep), pylab.NAN)
        return m
    
    def HeatMap(self, figure):
        lower_mat = self.MakeMatrix(self.lower_p)
        upper_mat = self.MakeMatrix(self.greater_p)
        
        nx = float(len(self.all_ind))
        ny = float(len(self.all_dep))
        x_indices = pylab.arange(nx) / nx  
        y_indices = pylab.arange(ny) / ny
        
        pylab.subplot(211, figure=figure)
        pylab.imshow(lower_mat, interpolation="nearest", figure=figure)
        pylab.title('Lower p-values', figure=figure)
        pylab.xticks(x_indices, self.all_ind, figure=figure)
        pylab.yticks(y_indices, self.all_dep, figure=figure)
        pylab.axis('scaled')
        
        pylab.subplot(212, figure=figure)
        pylab.imshow(upper_mat, interpolation="nearest", figure=figure)
        pylab.title('Upper p-values', figure=figure)
        pylab.xticks(x_indices, self.all_ind, figure=figure)
        pylab.yticks(y_indices, self.all_dep, figure=figure)
        pylab.axis('scaled')