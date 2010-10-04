class Bag(object):
    def __add__(self, other):
        result = self.copy()
        for item, count in other.iteritems():
            result._items[item] = result._items.get(item, 0) + count
        return result
    def __and__(self, other):
        result = Bag()
        for item, count in other.iteritems():
            new_count = min(self._items.get(item, 0), count)
            if new_count > 0:
                result._items[item] = new_count
        return result
    def __contains__(self, item):
        return item in self._items
    def __div__(self, integer):
        result = Bag(self)
        result.__idiv__(integer)
        return result
    def __getitem__(self, item):
        if (item not in self._items):
            self._items[item] = 0
        return self._items[item]            
    def __iadd__(self, other):
        self._items = self.__add__(other)._items
        return self
    def __iand__(self, other):
        self._items = self.__and__(other)._items
        return self
    def __idiv__(self, integer):
        for item in self.keys():
            self[item] /= integer
        return self
    def __init__(self, iterable=None):
        self._items = {}
        if iterable is not None:
            for item in iterable:
                self._items[item] = self._items.get(item, 0) + 1
    def __ior__(self, other):
        self._items = self.__or__(other)._items
        return self
    def __isub__(self, other):
        self._items = self.__sub__(other)._items
        return self
    def __iter__(self):
        for item, count in self.iteritems():
            for counter in xrange(count):
                yield item
    def __ixor__(self, other):
        self._items = self.__xor__(other)._items
        return self
    def __len__(self):
        return sum(self._items.itervalues())
    def __or__(self, other):
        result = self.copy()
        for item, count in other.iteritems():
            result._items[item] = max(result._items.get(item, 0), count)
        return result
    def to_list(self):
        result = []
        for item, count in self.iteritems():
            result += [item] * count
        return result
    def toset(self):
        return set(self.keys())
    def __repr__(self):
        return 'bag([%s])' % ', '.join([repr(item) for item in self.to_list()])
    def __setitem__(self, item, count):
        if not isinstance(count, int) or count < 0:
            raise ValueError
        if count > 0:
            self._items[item] = count
        elif item in self._items:
            del self._items[item]
    def __sub__(self, other):
        result = Bag()
        for item, count in self.iteritems():
            new_count = count - other._items.get(item, 0)
            if new_count > 0:
                result._items[item] = new_count
        return result
    def __xor__(self, other):
        result = self.copy()
        for item, count in other.iteritems():
            new_count = abs(result._items.get(item, 0) - count)
            if new_count > 0:
                result._items[item] = new_count
            elif item in result._item:
                del result._items[item]
        return result
    ### Comparative methods
    def issubset(self, other):
        for item, count in self.iteritems():
            if (count > other._items.get(item, 0)):
                return False
        return True
    def __lt__(self, other):
        return self.issubset(other) and (not other.issubset(self))
    def __le__(self, other):
        return self.issubset(other)
    def __gt__(self, other):
        return other.issubset(self) and (not self.issubset(other))
    def __ge__(self, other):
        return other.issubset(self)
    def __eq__(self, other):
        return self.issubset(other) and other.issubset(self)
    def __ne__(self, other):
        """
            Note that since we are dealing with multisets, there is no full ordering, which means
            that both a <= b and b <= a can be false.
        """
        return (not self.issubset(other)) or (not other.issubset(self))
    ### All methods
    def add(self, item, count=1):
        self._items[item] = self._items.get(item, 0) + count
    def discard(self, item):
        if (item in self._items.keys()):
            del self._items[item]
    def clear(self):
        self._items = {}
    def copy(self):
        result = Bag()
        result._items = self._items.copy()
        return result
    def difference(self, other):
        return self.__sub__(other)
    def difference_update(self, other):
        self._items = self.__sub__(other)._items
    def get(self, item, default=0):
        return self._items.get(item, default)
    def intersection(self, other):
        return self.__and__(other)
    def intersection_update(self, other):
        self.__iand__(other)
    def items(self):
        return self._items.items()
    def iteritems(self):
        return self._items.iteritems()
    def iterkeys(self):
        return self._items.iterkeys()
    def itervalues(self):
        return self._items.itervalues()
    def isempty(self):
        return (self._items == {})
    def keys(self):
        return self._items.keys()
    def pop(self):
        item = self._items.keys()[0]
        self._items[item] -= 1
        if self._items[item] == 0:
            del self._items[item]
        return item
    def remove(self, item, count=1):
        new_count = self._items[item] - count
        if new_count > 0:
            self._items[item] = new_count
        else:
            del self._items[item]
    def symmetric_difference(self, other):
        return self.__xor__(other)
    def symmetric_difference_update(self, other):
        self.__ixor__(other)
    def union(self, other):
        return self.__or__(other)
    def update(self, other):
        self.__ior__(other)
    def values(self):
        return self._items.values()
