import logging
import matcher
import Levenshtein
import re

from django.db.models import Q
from gibbs import models


class RegexApproxMatcher(matcher.Matcher):
    """A matcher that runs regular expressions directly on the database."""
            
    def _PrepareExpression(self, query):
        """Converts the query into a regular expression.
                
        Args:
            query: the string search query.
        
        Returns:
            A regular expression string.
        """
        if not query:
            return None
    
        # Escape regex special characters in the input (search for them
        # literally). Also, we allow '-', ',', '+' and digits in addition to spaces.
        query = re.escape(query.strip().lower())
        query = re.sub('(\\\?\s)+', '[-+\s\d,]+', query)
        # We allow leading and trailing junk.
        return '.*%s.*' % query

    def _FindNameMatches(self, query):
        """Override database search."""
        expression = self._PrepareExpression(query)
        if not expression:
            return []
        
        matches = models.CommonName.objects.filter(name__iregex=expression)
        return matches[:10*self._max_results]

class EditDistanceMatcher(RegexApproxMatcher):
    """Does near-exact matching and then uses edit distance to score."""
    
    def _GetScore(self, query, match):
        """Custom edit-distance based scoring."""
        str_query = str(query)
        str_candidate = str(match.key)
        dist = float(Levenshtein.distance(str_query, str_candidate))
        max_len = float(max(len(str_query), len(str_candidate)))
        return (max_len - dist) / max_len
    
    def _FindNameMatches(self, query):
        """Override database search."""
        qlen = len(query)
        if qlen < 5:
            matches = models.CommonName.objects.filter(name__icontains=query)
            return matches[:5*self._max_results]
        
        midpoint = qlen / 2
        head_re = self._PrepareExpression(query[:midpoint])
        tail_re = self._PrepareExpression(query[midpoint:])
        
        matches = models.CommonName.objects.filter(
            Q(name__iregex=head_re) | Q(name__iregex=tail_re))
        return matches[:5*self._max_results]
    

class CascadingMatcher(matcher.Matcher):
    """A matcher that tries multiple matching strategies."""
    
    def __init__(self, max_results=10, min_score=0.0):
        matcher.Matcher.__init__(self, max_results, min_score)
        self._exact_matcher = matcher.Matcher(max_results, min_score)
        self._re_matcher = RegexApproxMatcher(max_results, min_score)
        self._ed_matcher = EditDistanceMatcher(max_results, min_score)
    
    def _SortAndClip(self, matches):
        matches.sort(key=lambda m: m.score, reverse=True)
        return matches[:self._max_results]
    
    def Match(self, query):
        """Override base matching implementation."""  
        matches = self._exact_matcher.Match(query)
        
        match_set = set(m.key for m in matches)
        re_matches = self._re_matcher.Match(query)
        for m in re_matches:
            if m.key not in match_set:
                matches.append(m)
        
        if len(matches) >= self._max_results:
            return self._SortAndClip(matches)
        
        match_set = set(m.key for m in matches)
        ed_matches = self._ed_matcher.Match(query)
        for m in ed_matches:
            if m.key not in match_set:
                matches.append(m)
        
        return matches[:self._max_results]