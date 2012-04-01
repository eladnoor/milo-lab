#!/usr/bin/python

import numpy as np
import json

class GroupVector(list):
    """A vector of groups."""
    
    def __init__(self, groups_data, iterable=None):
        """Construct a vector.
        
        Args:
            groups_data: data about all the groups.
            iterable: data to load into the vector.
        """
        self.groups_data = groups_data
        
        if iterable is not None:
            self.extend(iterable)
        else:
            for _ in xrange(len(self.groups_data.all_group_names)):
                self.append(0)
    
    def __str__(self):
        """Return a sparse string representation of this group vector."""
        group_strs = []
        gv_flat = self.Flatten()
        for i, name in enumerate(self.groups_data.GetGroupNames()):
            if gv_flat[i]:
                group_strs.append('%s x %d' % (name, gv_flat[i]))
        return " | ".join(group_strs)
    
    def __iadd__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] += other[i]
        return self

    def __isub__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] -= other[i]
        return self
            
    def __add__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result[i] = self[i] + other[i]
        return result

    def __sub__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result[i] = self[i] - other[i]
        return result
    
    def __eq__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            if self[i] != other[i]:
                return False
        return True
    
    def __nonzero__(self):
        for i in xrange(len(self.groups_data.all_group_names)):
            if self[i] != 0:
                return True
        return False
    
    def __mul__(self, other):
        try:
            c = float(other)
            return GroupVector(self.groups_data, [x*c for x in self])
        except ValueError:
            raise ValueError("A GroupVector can only be multiplied by a scalar"
                             ", given " + str(other))
        
    def NetCharge(self):
        """Returns the net charge."""
        return int(np.dot(self, self.groups_data.all_group_charges))
    
    def Hydrogens(self):
        """Returns the number of protons."""
        return int(np.dot(self, self.groups_data.all_group_hydrogens))

    def Magnesiums(self):
        """Returns the number of Mg2+ ions."""
        return int(np.dot(self, self.groups_data.all_group_mgs))
    
    def RemoveEpsilonValues(self, epsilon=1e-10):
        for i in range(len(self)):
            if abs(self[i]) < epsilon:
                self[i] = 0
    
    def ToJSONString(self):
        return json.dumps(dict([(i, x) for (i, x) in enumerate(self) if x != 0]))
    
    @staticmethod
    def FromJSONString(groups_data, s):
        v = [0] * groups_data.Count()
        for i, x in json.loads(s).iteritems():
            v[int(i)] = x
        return GroupVector(groups_data, v)
    
    def Flatten(self):
        if not self.groups_data.transformed:
            return tuple(self)
        
        # map all pseudoisomeric group indices to Biochemical group indices (which are fewer)
        # use the names of each group and ignore the nH, z and nMg.
        biochemical_group_names = self.groups_data.GetGroupNames()
        biochemical_vector = [0] * len(biochemical_group_names)
        for i, x in enumerate(self):
            group_name = self.groups_data.all_groups[i].name
            new_index = biochemical_group_names.index(group_name)
            biochemical_vector[new_index] += x
        return tuple(biochemical_vector)        
    