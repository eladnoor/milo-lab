#!/usr/bin/python

class cartesian_product(object):
    
    def __init__(self, list_of_lists=[], unique_values=False):
        self._list_of_lists = list_of_lists
        self._unique_values = unique_values
        self._num_lists = len(self._list_of_lists)
        self._lengths = [len(x) for x in self._list_of_lists]
        self._total_size = 1
        for l in self._lengths:
            self._total_size *= l
        
        if (self._num_lists == 0 or self._total_size == 0):
            self._indices = None # this is like a flag indicating that the iteration is over
        else:
            self._indices = [0] * len(self._list_of_lists)
        
    def __iter__(self):
        return self
    
    def __len__(self):
        return self._total_size
    
    def get_current_item(self):
        if (self._indices == None):
            return []
        return [self._list_of_lists[i][self._indices[i]] for i in range(self._num_lists)]
    
    def increment(self):
        i = self._num_lists - 1
        self._indices[i] += 1
        
        while (i > 0 and self._indices[i] == self._lengths[i]):
            self._indices[i-1] += 1
            self._indices[i] = 0
            i -= 1

        if (self._indices[0] == self._lengths[0]):
            self._indices = None
    
    def next(self):
        if (self._unique_values):
            # if the current item contains repetitions, skip to the next one
            while (0 < len(set(self.get_current_item())) < self._num_lists):
                self.increment()

        if (self._indices == None):
            raise StopIteration()
        
        current_item = self.get_current_item()
        self.increment()
        
        return current_item

###############################################################
def test():
    for i in cartesian_product([[0,1],[0,2],[0,3]], unique_values=True):
        print i
  
if (__name__ == '__main__'):
    test()
