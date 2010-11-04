import logging
from gibbs import models
from util import topk


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
    """A class that matches a string against the database.
    
    The base implementation does exact matching.
    """
    
    def __init__(self, max_results=10, min_score=0.0):
        """Initializes the Matcher.
        
        Args:
            scorer: a MatchScorer object for scoring.
            max_results: the maximum number of matches to return.
            min_score: the minimum match score to return.
        """
        self._max_results = max_results
        self._min_score = min_score
    
    def _AcceptQuery(self, query):
        """Accept or rejec expression = self._PrepareExpression(query)
        results = models.CommonName.objects.filter(name__iregex=expression)t the query.
        
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
    
    def _FindNameMatches(self, query):
        """Find all the matches for this query.
        
        Args:
            query: the query to match.
            
        Returns:
            A list of CommonName objects matching the query.
        """
        try:
            name = models.CommonName.objects.get(name__iexact=query)
            return [name]
        except Exception, msg:
            return []
    
    def _GetScore(self, query, match):
        """Get the score for a query-match pair.
        
        Args:
            query: the query string.
            match: the Match object.
        
        Returns:
            A score between 0.0 and 1.0.
        """
        query_len = float(len(query))
        candidate_len = float(len(str(match.key)))
        return query_len / candidate_len
    
    def _ScoreMatches(self, query, matches):
        """Set the match scores for all matches.
        
        Args:
            query: the query string.
            matches: a list of match objects with uninitialized scores.
        """
        for m in matches:
            m.score = self._GetScore(query, m)
    
    def _FilterMatches(self, matches):
        """Filter the match list for min score.
        
        Args:
            matches: an unfiltered list of match objects.
        """ 
        return [m for m in matches if m.score >= self._min_score]
        
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
        name_matches = self._FindNameMatches(processed_query)
        
        matches = [Match(nm, None, 0.0) for nm in name_matches]
        for m in matches:
            compounds = m.key.compound_set.all()
            if compounds:
                m.value = compounds[0]
        self._ScoreMatches(processed_query, matches)
        self._FilterMatches(matches)
        
        matches.sort(key=lambda m: m.score, reverse=True)
        return matches
        