#!/usr/bin/python

import json


class EnzymeSet(set):
    """Represents a set of enzymes"""
         
    def __init__(self, name):
        set.__init__(self)
        
        self.name = name

    @staticmethod
    def FromDict(d):
        """Load from JSON data.
        
        Args:
            dict: a dictionary representing the EnzymeSet.
        """
        n = d.get('NAME', None)
        s = EnzymeSet(n)
        s.update(d.get('ECS', []))
        
        return s
    

class Pathway(object):
    
    def __init__(self, name, enzyme_sets):
        self.name = name
        self.enzyme_sets = enzyme_sets
        
    @staticmethod
    def FromDict(d):
        """Create a Pathway from a dictionary."""
        name = d.get('NAME', None)
        enzyme_sets = [EnzymeSet.FromDict(es) for es in d.get('ENZYME_SETS')]        
        return Pathway(name, enzyme_sets)


def LoadPathways(filename):
    f = open(filename)
    data = json.load(f)
    
    paths = [Pathway.FromDict(d) for d in data]
    f.close()
    return paths


if __name__ == '__main__':
    LoadPathways('../data/genomics/glycolysis_pathways.json')
