#!/usr/bin/python

import approximate_matcher
import reaction_parser
import unittest


class TestReactionParser(unittest.TestCase):
    """Tests for matcher.Match"""
    
    library = {'CO2': True, 'H2O': True, 'H+': True, 'H2': True, 'O2': True,
               'Q': True, 'Pi': True, 'NAD+': True, 'NADH': True, 'QH2': True,
               'GDP': True, 'GTP': True, 'Sodium': True, 'hydrogen chloride': True,
               'Sodium Chloride': True, 'Acetyl-CoA': True, 'CoA-SH': True}
    
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
            self.assertNotEqual(None, parsed, msg=reaction_str)
            
            for compound in parsed.reactants + parsed.products:
                self.assertNotEqual(None, compound.parsed_name)
                self.assertNotEqual(None, compound.parsed_coeff)
                self.assertTrue(compound.parsed_coeff > 0)
                self.assertTrue(len(compound.matches) > 0)
                
if __name__ == '__main__':
    unittest.main()