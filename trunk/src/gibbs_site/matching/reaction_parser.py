import re
import logging
import pyparsing


def _parsedCompound(c_list):
    """Always put a stoichiometric coefficient with a compound."""
    if len(c_list) == 2:
        return c_list[0], c_list[1]
    return 1, c_list[0]


def _makeparser():
    """Builds a pyparsing-based recursive descent parser for chemical reactions."""
    coeff = pyparsing.Word(pyparsing.nums).setParseAction(lambda t:int(t[0]))
    optional_coeff = pyparsing.Optional(coeff)
    compound_separator = pyparsing.Literal('+').suppress()
    
    compound_name_component = pyparsing.Word(pyparsing.alphanums, pyparsing.alphanums + "-+,()'")
    compound_name = pyparsing.Forward()
    compound_name << (compound_name_component + pyparsing.ZeroOrMore(compound_name_component))
    compound_name.setParseAction(lambda s: ' '.join(s))
    
    compound_with_coeff = pyparsing.Forward()
    compound_with_coeff << ((optional_coeff + compound_name) | compound_name)
    compound_with_coeff.setParseAction(_parsedCompound)
    compound_with_coeff.setResultsName("compound")
    
    compound_with_separator = pyparsing.Forward()
    compound_with_separator << (compound_with_coeff + compound_separator)
    
    reaction_side = pyparsing.Forward()
    reaction_side << (pyparsing.ZeroOrMore(compound_with_separator) +
                      compound_with_coeff)
    reaction_side.setParseAction(lambda l: [l])
    reaction_side.setResultsName("reaction_side")
    
    side_separators = [pyparsing.Literal(s) for s in ("=", "->", "=>", "<=>")]
    side_separator = pyparsing.Or(side_separators).suppress()
    
    reaction = pyparsing.Forward()
    reaction << (reaction_side + side_separator + reaction_side)
    return reaction


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
    
    REACTION_PATTERN = r'.*\s+(=>|<=>|=|->|<->)\s+.*'
    REACTION_MATCHER = re.compile(REACTION_PATTERN)
    
    def __init__(self, compound_matcher):
        """Initialize the reaction parser.
        
        Args:
            compound_matcher: a matcher.Matcher object that matches
                              individual compounds.
        """
        self._matcher = compound_matcher
        self._parser = _makeparser()
    
    def ShouldParseAsReaction(self, query):
        """Returns True if this query is likely to be a reaction query.
        
        Args:
            query: the query string.
        """
        m = self.REACTION_MATCHER.match(query.strip())
        return m is not None
    
    def _MakeReactionCompoundMatch(self, coeff, name):
        compound_matches = self._matcher.Match(name)
        return ReactionCompoundMatch(name, coeff, compound_matches)
    
    def _MakeReactionSide(self, parsed_side):
        side = []
        for coeff, name in parsed_side:
            side.append(self._MakeReactionCompoundMatch(coeff, name))
        return side
    
    def ParseReactionQuery(self, query):
        """Parse the query as a reaction.
        
        Args:
            query: the query string.
        
        Returns:
            An initialized ParsedReaction object, or None if parsing failed.
        """
        results = None
        try:
            results = self._parser.parseString(query)
        except Exception, e:
            logging.error(e)
            return None
        
        reactants, products = results
        
        reactants = self._MakeReactionSide(reactants)
        products = self._MakeReactionSide(products)
        
        if not reactants:
            logging.warn('Failed to parse reactants.')
            return None
        
        if not products:
            logging.warn('Failed to parse products.')
            return None
        
        return ParsedReaction(reactants, products)
        
        

        