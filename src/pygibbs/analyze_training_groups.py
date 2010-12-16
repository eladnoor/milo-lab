#!/usr/bin/python

import logging
import os
import pylab

from pygibbs import group_decomposition
from pygibbs import groups_data
from pygibbs import pseudoisomers_data
from pygibbs import templates
from toolbox import util

class GroupMap(object):
    
    def __init__(self, map):
        self.map = map
        
        self.counts = dict((g, len(s)) for g,s in self.map.iteritems())

    def GetGroupsByNumExamples(self, num_examples):
        matches = [(g,pisomers) for g,pisomers in self.map.iteritems()
                   if len(pisomers) == num_examples]
        return matches
    
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
    
    def __init__(self, groups_data):
        self.groups_data = groups_data
        self.decompositions = {}
        
    def AddDecomposition(self, pisomer, decomposition):
        self.decompositions[pisomer] = decomposition
    
    def GetGroupMap(self):
        map = {}
        for pisomer, decomposition in self.decompositions.iteritems():
            for group, node_sets in decomposition.SparseRepresentation().iteritems():
                map.setdefault(group, set()).add(pisomer)
        
        return GroupMap(map)
                

def main():
    gdata = groups_data.GroupsData.FromGroupsFile(
        '../data/thermodynamics/groups_species.csv')
    decomposer = group_decomposition.GroupDecomposer(gdata)
    pdata = pseudoisomers_data.PseudoisomersData.FromFile(
        '../data/thermodynamics/dG0.csv')
    
    dstats = DecompositionStats(groups_data)

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
        
    template_data = {'rare_groups_by_count': rare_groups,
                     'most_common_groups': mcg_dicts,
                     'histo_image_name': image_name}
    templates.render_to_file('analyze_training_groups.html',
                             template_data,
                             '../res/analyze_training_groups.html')
    
    
    
if __name__ == '__main__':
    main()
    logging.info('Done')