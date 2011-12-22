#!/usr/bin/python
    
import numpy as np


class Bounds(object):
    """Contains upper and lower bounds for various keys."""
    
    def __init__(self,
                 lower_bounds=None,
                 upper_bounds=None,
                 default_lb=None,
                 default_ub=None):
        """Initialize the bounds object.
        
        Args:
            lower_bounds: a dictionary mapping strings to float lower bounds.
            upper_bounds: a dictionary mapping strings to float upper bounds.
            default_lb: the default lower bound to return.
            default_lb: the default upper bound to return.
        """
        self.lower_bounds = lower_bounds or {}
        self.upper_bounds = upper_bounds or {}
        self.default_lb = default_lb
        self.default_ub = default_ub
        
    
    def GetLowerBound(self, key):
        """Get the lower bound for this key.
        
        Args:
            key: a string key.
        """
        if self.lower_bounds and key in self.lower_bounds:
            return self.lower_bounds[key]
        return self.default_lb
    
    def GetUpperBound(self, key):
        """Get the upper bound for this key.
        
        Args:
            key: a string key.
        """
        if self.upper_bounds and key in self.upper_bounds:
            return self.upper_bounds[key]
        return self.default_ub
    
    def GetBounds(self, keys):
        """Get the bounds for a set of keys in order.
        
        Args:
            keys: an iterable of keys.
        
        Returns:
            A two-tuple (lower_bounds, upper_bounds)
        """
        lower_bounds = np.array([self.GetLowerBound(key) for key in keys])
        upper_bounds = np.array([self.GetUpperBound(key) for key in keys])
        return lower_bounds, upper_bounds

    def GetLnBounds(self, keys):
        """Get the bounds for a set of keys in order.
        
        Args:
            keys: an iterable of keys.
        
        Returns:
            A two-tuple (lower_bounds, upper_bounds)
        """
        lb, ub = self.GetBounds(keys)
        return np.log(lb), np.log(ub)
    