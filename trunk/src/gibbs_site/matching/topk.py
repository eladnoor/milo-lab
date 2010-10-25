"""A data structure that keeps only the top K elements it sees."""

import heapq

class TopK(object):
    """Keeps the K top items."""
    
    def __init__(self, max):
        """Construction.
        
        Args:
            max: the maximum number of items to keep (k).
        """
        self._heap = []
        self.max = max
    
    def _Smallest(self):
        """Get the smallest item."""
        heapq.nsmallest(1, self._heap)[0]
        
    def GetSorted(self, key=None):
        """Return top items as a sorted list using the key function."""
        return heapq.nlargest(self.max, self._heap, key=key)

    def MaybeAdd(self, elt):
        """Potentially add elt if it's big enough.
        
        Args:
            elt: the element to add. must implement comparison.
        """
        if len(self._heap) < self.max:
            heapq.heappush(self._heap, elt)
            return
        
        if self._Smallest() < elt:
            heapq.heappop(self._heap)
            heapq.heappush(self._heap, elt)