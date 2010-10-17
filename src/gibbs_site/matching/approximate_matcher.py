import matcher
import Levenshtein
import re


class EditDistanceMatcher(matcher.Matcher):
    """A linear-scanning matcher that does approximate string matching.
    
    Uses an edit-distance implementation to do approximate matching.
    """
    
    def _MatchSingle(self, query, candidate):
        """Override the single-candidate matching implementation."""
        str_query = str(query)
        str_candidate = str(candidate)
        dist = float(Levenshtein.distance(str_query, str_candidate))
        max_len = float(max(len(str_query), len(str_candidate)))
        return (max_len - dist) / max_len


class RegexApproxMatcher(matcher.Matcher):
    """A linear-scanning matcher that does approximate string matching.
    
    Uses regular-expression approach to do approximate matching.
    """
    
    _COMPILED_EXPRESSION_CACHE = {}
        
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
    
    def _GetCompiledExpression(self, query):
        s = query.strip()
        if not s:
            return None
        
        if s in self._COMPILED_EXPRESSION_CACHE:
            return self._COMPILED_EXPRESSION_CACHE[s]
        
        re_expression = self._PrepareExpression(query)
        if not re_expression:
            return None
        
        compiled_re = re.compile(re_expression)
        self._COMPILED_EXPRESSION_CACHE[s] = compiled_re
        return compiled_re
        
    def _MatchSingle(self, query, candidate):
        """Override the single-candidate matching implementation."""        
        compiled_re = self._GetCompiledExpression(query)
        if not compiled_re:
            return 0.0
        
        m = compiled_re.search(candidate)
        if not m:
            return 0.0
        
        query_len = float(len(query))
        candidate_len = float(len(candidate))
        return query_len / candidate_len
    
    
class BackfillingRegexApproxMatcher(matcher.Matcher):
    """A matcher that uses Regex by default and falls back to
       Levenshtein distance when there are not enough matches.
    """
    
    def __init__(self, library, max_results=10, min_score=0.0):
        matcher.Matcher.__init__(self, library, max_results, min_score)
        self._re_matcher = RegexApproxMatcher(library, max_results, min_score)
        self._edit_dist_matcher = EditDistanceMatcher(library, max_results, min_score)
    
    def Match(self, query):
        """Override base matching implementation."""
        matches = self._re_matcher.Match(query)
        if len(matches) == self._max_results:
            return matches
        
        match_set = set(m.key for m in matches)
        ed_matches = self._edit_dist_matcher.Match(query)
        for m in ed_matches:
            if m.key not in match_set:
                matches.append(m)
        
        return matches[:self._max_results]