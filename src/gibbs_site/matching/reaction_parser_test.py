#!/usr/bin/python

import approximate_matcher
import reaction_parser
import unittest


class TestReactionParser(unittest.TestCase):
    """Tests for matcher.Match"""
    
    library = ('CO2', 'H2O', 'H+', 'H2', 'O2', 'Q', 'Pi',
               'NAD+', 'NADH', 'QH2', 'GDP', 'GTP',  
               'Sodium', 'hydrogen chloride', 'Sodium Chloride',
               'Acetyl-CoA', 'CoA-SH')
    
    parsable_reactions = (
        'H2O => 2 H2 + O2',
        '2 Sodium + 2 Hydrogen chloride <=> 2 Sodium chloride + H2',
        'Acetyl-CoA + 3 NAD+ + Q + GDP + Pi + 2 H2O <=> '
        'CoA-SH + 3 NADH + 3 H+ + QH2 + GTP + 2 CO2'
        )
    
    def setUp(self):
        self._matcher = approximate_matcher.RegexApproxMatcher(
            self.library, max_results=1, min_score=0.4)
        self._parser = reaction_parser.ReactionParser(self._matcher)
    
    def testParsing(self):
        for reaction_str in self.parsable_reactions:
            self.assertTrue(self._parser.ShouldParseAsReaction(reaction_str))
            parsed = self._parser.ParseReactionQuery(reaction_str)
            self.assertNotEqual(None, parsed)
            
            for compound in parsed.reactants + parsed.products:
                self.assertNotEqual(None, compound.parsed_name)
                self.assertNotEqual(None, compound.parsed_coeff)
                self.assertTrue(compound.parsed_coeff > 0)
                self.assertTrue(len(compound.matches) > 0)
                
if __name__ == '__main__':
    unittest.main()