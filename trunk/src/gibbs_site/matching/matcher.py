class Error(Exception):
    pass


class IllegalQueryError(Error):
    pass


class Match(object):
    """An object containing a string match and it's score."""
    
    def __init__(self, key, value, score):
        """Initialize a Match.
        
        Args:
            key: the object that matched.
            value: the value of the match (the object pointed to by the key).
            score: between 0.0 and 1.0, higher is better.
        """
        self.key = key
        self.value = value
        self.score = score
    
    def __eq__(self, other):
        """Equality checking between matches, used for testing."""
        return (self.key == other.key and
                self.value == other.value and
                self.score == other.score)
    
    def __str__(self):
        """Get as a string for debugging/printability."""
        return '<matcher.Match> value=%s, score=%f' % (self.value,
                                                       self.score)


class Matcher(object):
    """A class that matches a string against a library of strings.
    
    The base implementation does exact matching.
    """
    
    def __init__(self, library, max_results=10, min_score=0.0):
        """Initializes the Matcher.
        
        Args:
            library: a dictionary mapping names (to match against) to objects (to return).
            max_results: the maximum number of matches to return.
            min_score: the minimum match score to return.
        """
        self._library = library
        self._max_results = max_results
        self._min_score = min_score
    
    def _AcceptQuery(self, query):
        """Accept or reject the query.
        
        Returns:
            True if the query is accepted.
        """
        if query.strip():
            return True
        return False
    
    def _PreprocessQuery(self, query):
        """Perform pre-search query manipulation.
        
        Default implementation simply strips leading/trailing whitespace
        and makes the query lowercase.
        
        Args:
            query: the string query.
            
        Returns:
            The pre-processed query as a string.
        """
        return query.strip().lower()
    
    def _PrepocessCandidate(self, candidate):
        """Perform pre-match candidate manipulation.
        
        Default implementation converts to a lower-case string.
        
        Args:
            candidate: the candidate object (convertible to a string).
        
        Returns:
            The pre-processed candidate as a string.
        """
        return str(candidate).strip().lower()
    
    def _MatchSingle(self, query, candidate):
        """Checks a single query, candidate pair for a match.
        
        For linear-scanning matchers, it suffices to override this method.
        Base implementation is exact matching.
        
        Args:
            query: the query string.
            candidate: the candidate match string.
        
        Returns:
            A score between 0.0 and 1.0 for the match between query and
            candidate. Higher is better.
        """
        if query.strip() == candidate.strip():
            return 1.0
        return 0.0
    
    def Match(self, query):
        """Find matches for the query in the library.
        
        Args:
            query: the string query.
        
        Returns:
            A sorted list of Match objects or None if
            the query could not be parsed.
        """
        if not self._AcceptQuery(query):
            raise IllegalQueryError('%s is not a valid query' % query)
        
        processed_query = self._PreprocessQuery(query)
        
        matches = []
        for candidate in self._library.iterkeys():
            score = self._MatchSingle(processed_query,
                                      self._PrepocessCandidate(candidate))
            if score > self._min_score:
                matches.append(Match(candidate, self._library[candidate], score))
        
        # Sort descending by score and return the top k.
        matches.sort(key=lambda x: x.score, reverse=True)
        return matches[:self._max_results]