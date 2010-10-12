#!/usr/bin/python

import matcher
import unittest


class TestMatch(unittest.TestCase):
    """Tests for matcher.Match"""
    
    def testConstruction(self):
        m = matcher.Match('a', 0.1)
        self.assertEqual('a', m.value)
        self.assertEqual(0.1, m.score)


class FirstLastCharacterMatcher(matcher.Matcher):
    """A test matcher.
    
    Considers a match at the first character perfect and a match
    at the last character slightly worse.
    """
    
    def _MatchSingle(self, query, candidate):
        """Override the single-candidate matching implementation."""
        if query[0] == candidate[0]:
            return 1.0
        if query[-1] == candidate[-1]:
            return 0.75
        return 0.0


class TestMatcher(unittest.TestCase):
    """Tests for matcher.Matcher."""
    
    library = ('avi',
               'ron',
               'eran',
               'elad',
               'niv',
               'shira',
               'lior',
               'libat')
    
    def _CheckIsSortedByScore(self, results):
        prev_score = 10.0 # Scores are all <= 1.0
        for match in results:
            self.assertTrue(match.score <= prev_score)
            prev_score = match.score
    
    def testBaseMatcher(self):
        m = matcher.Matcher(self.library, 10)
        
        for name in self.library:
            results = m.Match(name)
            self.assertEqual(1, len(results))
            self.assertEqual(matcher.Match(name, 1.0), results[0])
    
    def testSortingMatcher(self):
        m = FirstLastCharacterMatcher(self.library, 10)
        
        results = m.Match('ron')
        self.assertEqual(2, len(results))
        self._CheckIsSortedByScore(results)
        
        results = m.Match('ligand')
        self.assertEqual(3, len(results))
        self._CheckIsSortedByScore(results)
        
    def testMaxResults(self):
        m = FirstLastCharacterMatcher(self.library, 2)
        
        results = m.Match('ligand')
        self.assertEqual(2, len(results))
        self._CheckIsSortedByScore(results)
    
    def testMinScore(self):
        m = FirstLastCharacterMatcher(self.library,
                                      max_results=10,
                                      min_score=0.8)
        
        results = m.Match('ligand')
        self.assertEqual(2, len(results))
        self._CheckIsSortedByScore(results)
    
    
if __name__ == '__main__':
    unittest.main()