import re
import logging
import matcher


class ReactionCompoundMatch(object):
    """A match of a compound in a reaction.
    
    Contains the parsed name (what the user entered), the parsed
    stoichiometric coefficient, and a list of matcher.Match objects
    of the potential matches.
    """
    def __init__(self, parsed_name, parsed_coeff, matches):
        self.parsed_name = parsed_name
        self.parsed_coeff = parsed_coeff
        self.matches = matches
    
    def __str__(self):
        return '%d %s, matches: %s' % (self.parsed_coeff,
                                       self.parsed_name,
                                       ', '.join(self.matches))
    
    def ParsedDataEqual(self, other):
        """Checks if the parsed data (name, coefficient) are equal
           for two ReactionCompoundMathes.
        """
        return (self.parsed_coeff == other.parsed_coeff and
                self.parsed_name == other.parsed_name)

    
class ParsedReaction(object):
    """A reaction parsed from a query."""
    
    def __init__(self, reactants=None, products=None):
        """Initialize the ParsedReaction object.
        
        Args:
            reactants: a list of ReactionCompoundMatches for the reactants.
            products: a list of ReactionCompoundMatches for the products.
        """
        self.reactants = reactants or []
        self.products = products or []
    
    def GetBestMatch(self):
        """Returns the product and reactant lists for the best match.
        
        Each list is of 2 tuples (coeff, kegg_id).
        """
        reactants = []
        for c in self.reactants:
            if not c.matches:
                return None, None
            reactants.append((c.parsed_coeff, c.matches[0].value.kegg_id))
        
        products = []
        for c in self.products:
            if not c.matches:
                return None, None
            products.append((c.parsed_coeff, c.matches[0].value.kegg_id))
        
        return reactants, products
    

class ReactionParser(object):
    """Parses reaction queries from users."""
    
    REACTION_PATTERN = r'.*\s+(=>|<=>|=)\s+.*'
    REACTION_MATCHER = re.compile(REACTION_PATTERN)
    
    REACTION_SPLITTER = re.compile(r'\s*(?:=>|<=>|=)\s*')
    REACTION_SIDE_SPLITTER = re.compile(r'(?:\s+\+\s+)')
    
    COMPOUND_COEFF_MATCHER = re.compile(r'^(\d+\s+)?\s*(.*)$')
    
    def __init__(self, compound_matcher):
        """Initialize the reaction parser.
        
        Args:
            compound_matcher: a matcher.Matcher object that matches
                              individual compounds.
        """
        self._matcher = compound_matcher
    
    def ShouldParseAsReaction(self, query):
        """Returns True if this query is likely to be a reaction query.
        
        Args:
            query: the query string.
        """
        m = self.REACTION_MATCHER.match(query.strip())
        return m is not None
    
    def _ParseCompoundAndCoefficient(self, compound_str):
        """Parses the string containing a compound and coefficient.
        
        Args:
            compound_str: the string containing a compound and its
                          stoichiometric coefficient.
        
        Returns:
            A ReactionCompoundMatch or None if parsing failed.
        """
        m = self.COMPOUND_COEFF_MATCHER.match(compound_str)
        if not m:
            return None
        
        coeff, name = m.groups()
        coeff = coeff or '1'
        coeff = int(coeff.strip())
            
        compound_matches = self._matcher.Match(name)
        return ReactionCompoundMatch(name, coeff, compound_matches)
    
    def _ParseReactionSide(self, side_list):
        """Parses a side of a reaction.
        
        Args:
            side_list: the list of strings of compounds with stoichiometric
                       coefficients.
        
        Returns:
            A list of ReactionCompoundMatches or None if any one failed to parse.
        """
        side_compounds = []
        for compound_str in side_list:
            parsed_compound = self._ParseCompoundAndCoefficient(compound_str)
            if not parsed_compound:
                return None
            side_compounds.append(parsed_compound)
        
        return side_compounds
        
    def ParseReactionQuery(self, query):
        """Parse the query as a reaction.
        
        Args:
            query: the query string.
        
        Returns:
            An initialized ParsedReaction object, or None if parsing failed.
        """
        split_query = self.REACTION_SPLITTER.split(query)
        if len(split_query) != 2:
            logging.warn('Failed to parse reaction query: "%s"', query)
            return None
        
        reactants_str, products_str = split_query
        reactants_str = reactants_str.strip()
        products_str = products_str.strip()
        
        reactants_with_coeffs = self.REACTION_SIDE_SPLITTER.split(reactants_str)
        products_with_coeffs = self.REACTION_SIDE_SPLITTER.split(products_str)
        
        reactants = self._ParseReactionSide(reactants_with_coeffs)
        products = self._ParseReactionSide(products_with_coeffs)
        
        if not reactants:
            logging.warn('Failed to parse reactants.')
            return None
        
        if not products:
            logging.warn('Failed to parse products.')
            return None
        
        return ParsedReaction(reactants, products)
        
        

        