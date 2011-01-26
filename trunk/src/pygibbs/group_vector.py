#!/usr/bin/python

import pylab

class GroupVector(list):
    """A vector of groups."""
    
    def __init__(self, groups_data, iterable=None):
        """Construct a vector.
        
        Args:
            groups_data: data about all the groups.
            iterable: data to load into the vector.
        """
        self.groups_data = groups_data
        
        if iterable:
            self.extend(iterable)
    
    def __str__(self):
        """Return a string representation of this group vector."""
        group_strs = []
        
        for i, name in enumerate(self.groups_data.all_group_names):
            if self[i]:
                group_strs.append('%s x %d' % (name, self[i]))
        return " | ".join(group_strs)
    
    def __iadd__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] += other[i]

    def __isub__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] -= other[i]
            
    def __add__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result.append(self[i] + other[i])
        return result

    def __sub__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result.append(self[i] - other[i])
        return result
    
    def __eq__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            if self[i] != other[i]:
                return False
        return True 
    
    def NetCharge(self):
        """Returns the net charge."""
        return int(pylab.dot(self, self.groups_data.all_group_charges))
    
    def Hydrogens(self):
        """Returns the number of protons."""
        return int(pylab.dot(self, self.groups_data.all_group_hydrogens))

    def Magnesiums(self):
        """Returns the number of Mg2+ ions."""
        return int(pylab.dot(self, self.groups_data.all_group_mgs))