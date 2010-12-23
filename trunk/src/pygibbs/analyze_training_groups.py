#!/usr/bin/python

import logging
import os
import pylab

from pygibbs import group_decomposition
from pygibbs import groups_data
from pygibbs import pseudoisomers_data
from pygibbs import templates
from toolbox import util
from gibbs_site.util import topk


class PairwiseGroupMap(object):
    
    def __init__(self, decompositions, group_map, pairwise_map):
        self.map = pairwise_map
        self.group_map = group_map
        self.decompositions = decompositions
        
    def GetMostFrequentCooccurences(self, num=20):
        tk = topk.TopK(num)
        for groups, pisomers in self.map.iteritems():
            ga, gb = groups
            tk.MaybeAdd((len(pisomers), ga, gb))
        
        return [(ga, gb, count) for count, ga, gb in
                tk.GetSorted(key=lambda s:s[0])]
    
    def _GetRatio(self, pisomer, groupa, groupb):
        sparse = self.decompositions[pisomer].SparseRepresentation()
        return float(len(sparse[groupa])) / float(len(sparse[groupb]))
    
    def GetLinearlyDependentPairs(self):
        s = set()
        for groups, pisomers in self.map.iteritems():            
            ga, gb = groups
            ratio = None
            ratio_matches = True
            for pisomer in pisomers:
                current_ratio = self._GetRatio(pisomer, ga, gb)
                
                if not ratio:
                    ratio = current_ratio
                    continue
                
                if ratio != current_ratio:
                    ratio_matches = False
                    break
                
            if not ratio_matches:
                continue
            
            num_together = len(pisomers)
            counts = self.group_map.GetCounts()
            acount, bcount = counts[ga], counts[gb]
            if not num_together == acount or not num_together == bcount:
                continue
            
            s.add((ga, gb, ratio, frozenset(pisomers)))
        
        return s

class GroupMap(object):
    
    def __init__(self, map):
        self.map = map
        
        self.counts = dict((g, len(s)) for g,s in self.map.iteritems())

    def GetCounts(self):
        return self.counts
    
    def AverageNumExamples(self):
        return pylab.mean(self.counts.values())
    
    def MedianNumExamples(self):
        return pylab.median(self.counts.values())
    
    def StdDevNumExamples(self):
        return pylab.std(self.counts.values())
    
    def NumExamplesRange(self):
        cts = self.counts.values()
        return min(cts), max(cts)
    
    average_examples = property(AverageNumExamples)
    median_examples = property(MedianNumExamples)
    std_dev_examples = property(StdDevNumExamples)
    examples_range = property(NumExamplesRange)

    def GetGroupsByNumExamples(self, num_examples):
        matches = [(g,pisomers) for g,pisomers in self.map.iteritems()
                   if len(pisomers) == num_examples]
        return matches
    
    zeroes = property(lambda self: self.GetGroupsByNumExamples(0))
    ones = property(lambda self: self.GetGroupsByNumExamples(1))
    
    def GetMostCommonGroups(self, num_groups=10):
        counts_per_group = [(c, g) for g,c in self.counts.iteritems()]
        counts_per_group.sort(reverse=True)
        return [(g,c) for c,g in counts_per_group[:num_groups]]
    
    def PlotHistogram(self, filename):
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            util._mkdir(dirname)
        
        fig = pylab.figure()
        
        c = list(self.counts.values())
        pylab.title('Count Per Group Histogram')
        pylab.xlabel('Count')
        pylab.ylabel('Number of Groups')
        pylab.hist(c, pylab.arange(0, max(c), 1))
        fig.savefig(filename, format='png')
    

class DecompositionStats(object):
    
    def __init__(self, gdata):
        self.groups_data = gdata
        self.decompositions = {}
        
    def AddDecomposition(self, pisomer, decomposition):
        self.decompositions[pisomer] = decomposition
    
    def GetGroupMap(self):
        map = dict((group, set()) for group in self.groups_data.all_groups)
        for pisomer, decomposition in self.decompositions.iteritems():
            for group, node_sets in decomposition.SparseRepresentation().iteritems():
                map.setdefault(group, set()).add(pisomer)
        
        return GroupMap(map)
    
    def GetPairwiseGroupMap(self):
        group_cooccurrences = {}
        for pisomer, decomposition in self.decompositions.iteritems():
            sparse = decomposition.SparseRepresentation()
            for ga in sparse:
                for gb in sparse:
                    if ga != gb:
                        key = (ga, gb)
                        if (gb, ga) in group_cooccurrences:
                            key = (gb, ga)
                        group_cooccurrences.setdefault(key, set()).add(pisomer)
        return PairwiseGroupMap(self.decompositions,
                                self.GetGroupMap(),
                                group_cooccurrences)
                
                

def main():
    gdata = groups_data.GroupsData.FromGroupsFile(
        '../data/thermodynamics/groups_species.csv')
    decomposer = group_decomposition.GroupDecomposer(gdata)
    pdata = pseudoisomers_data.PseudoisomersData.FromFile(
        '../data/thermodynamics/dG0.csv')
    
    dstats = DecompositionStats(gdata)

    for pisomer in pdata:
        if not pisomer.Train():
            continue
        
        mol = pisomer.Mol()
        decomposition = None
        if mol:
            decomposition = decomposer.Decompose(pisomer.Mol())
            dstats.AddDecomposition(pisomer, decomposition)
        else:
            logging.warning('Cannot get a Mol for pseudoisomer %s' % pisomer)
            continue
        
    map = dstats.GetGroupMap()
    rare_groups = []
    for i in range(10):
        groups = map.GetGroupsByNumExamples(i)
        l = [{'group': g, 'pseudoisomers': ps} for g,ps in groups]
        rare_groups.append({'count': i, 'groups': l})
    
    image_name = 'images/groups_histo.png'
    map.PlotHistogram('../res/' + image_name)
    
    most_common_groups = map.GetMostCommonGroups()
    mcg_dicts = [{'count': c, 'group': g} for g,c in most_common_groups]
    
    pairwise_map = dstats.GetPairwiseGroupMap()
    frequent_pairs = pairwise_map.GetMostFrequentCooccurences()
    freq_pairs_dicts = [{'count': c, 'groups': [ga, gb]} for ga, gb, c in frequent_pairs]
    ld_pairs = pairwise_map.GetLinearlyDependentPairs()
    ld_pairs_dicts = [{'ratio': r, 'pisomers': p, 'count': len(p), 'groups': [ga, gb]}
                      for ga, gb, r, p in ld_pairs]
    
    template_data = {'rare_groups_by_count': rare_groups,
                     'most_common_groups': mcg_dicts,
                     'histo_image_name': image_name,
                     'frequent_pairs': freq_pairs_dicts,
                     'ld_pairs': ld_pairs_dicts,
                     'groups_data': gdata,
                     'group_map': map}
    templates.render_to_file('analyze_training_groups.html',
                             template_data,
                             '../res/analyze_training_groups.html')
    
    
    
if __name__ == '__main__':
    main()
    logging.info('Done')