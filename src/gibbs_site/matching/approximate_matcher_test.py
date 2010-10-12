#!/usr/bin/python

import approximate_matcher
import unittest

class TestMatcher(unittest.TestCase):
    """Tests for matcher.Matcher."""
    
    library = ('avi',
               'ari',
               'ariel',
               'glucose',
               'glucosamine',
               'alanine',
               'phenylalanine',
               'l-glucosamine')
    
    def _CheckIsSortedByScore(self, results):
        prev_score = 10.0 # Scores are all <= 1.0
        for match in results:
            self.assertTrue(match.score <= prev_score)
            prev_score = match.score
    
    def testEditDistanceMatcher(self):
        m = approximate_matcher.EditDistanceMatcher(self.library,
                                                    max_results=10,
                                                    min_score=0.3)
    
        results = m.Match('avi')
        self.assertEqual(3, len(results))
        self._CheckIsSortedByScore(results)
        
        results = m.Match('d-glucose')
        self.assertEqual(3, len(results))
        self._CheckIsSortedByScore(results)
        
        results = m.Match('acetyl-coa')
        self.assertEqual(0, len(results))

    def testPrepareExpression(self):
        examples = (('  teSt    tEsT ', '.*test[-+\s]test.*'),
                    ('gluco', '.*gluco.*'),
                    ('D Fructo', '.*d[-+\s]fructo.*'))
        for query, expression in examples:
            self.assertEqual(expression,
                             approximate_matcher.RegexApproxMatcher._PrepareExpression(query))

    def testRegexApproxMatcher(self):
        m = approximate_matcher.RegexApproxMatcher(self.library,
                                                   max_results=10,
                                                   min_score=0.3)
        
        # Has an exact match.
        results = m.Match('avi')
        self.assertEqual(1, len(results))
        self._CheckIsSortedByScore(results)
        
        # Has no matches because the 'd=' doesn't appear in the library. 
        results = m.Match('d-glucose')
        self.assertEqual(0, len(results))
        
        # Nothing similar in the library.
        results = m.Match('acetyl-coa')
        self.assertEqual(0, len(results))
        
        # Has several substring matches.
        results = m.Match('gluco')
        self.assertEqual(3, len(results))
        self._CheckIsSortedByScore(results)

        # Has several substring matches.
        results = m.Match('alan')
        self.assertEqual(2, len(results))
        self._CheckIsSortedByScore(results)  

if __name__ == '__main__':
    unittest.main()