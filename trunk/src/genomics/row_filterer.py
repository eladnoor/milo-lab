#!/usr/bin/python


class RowFilterer(object):
    """Filters rows (dictionaries)."""
    def __init__(self, filter_cols, filter_vals):
        """Initialization.
        
        Args:
            filter_cols: the columns (dictionary keys) to filter.
            filter_vals: the column values to keep.
        """
        self.cols = filter_cols
        self.vals = filter_vals
        
    def Keep(self, row):
        """Returns true if the row should be kept.
        
        Args:
            row: a dictionary of a row (most likely from a DictReader).
        
        Returns:
            True if the row should be kept.
        """
        # Apply column filters
        apply_filter = lambda x, y: x in row and row[x] == y
        applied_filters = map(apply_filter, self.cols, self.vals)
        or_reduce = lambda x, y: x or y
        passed_filter = reduce(or_reduce, applied_filters, False)
        return passed_filter