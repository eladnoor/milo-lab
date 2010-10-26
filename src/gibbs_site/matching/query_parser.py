import logging
import pyparsing
import re

def _parsedCompound(c_list):
    """Always put a stoichiometric coefficient with a compound."""
    if len(c_list) == 2:
        return c_list[0], c_list[1]
    return 1, c_list[0]


def _MakeReactionParser():
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


class ParsedReactionQuery(object):
    """A parsed reaction query."""
    
    def __init__(self, reactants=None, products=None):
        """Initialize the ParsedReaction object.
        
        Args:
            reactants: a list of ReactionCompoundMatches for the reactants.
            products: a list of ReactionCompoundMatches for the products.
        """
        self.reactants = reactants or []
        self.products = products or []
    

class QueryParser(object):
    
    REACTION_PATTERN = r'.*(=>|<=>|=|->|<->).*'
    REACTION_MATCHER = re.compile(REACTION_PATTERN)
    
    def __init__(self):
        """Initialize the parser."""
        self._rparser = _MakeReactionParser()
        
    def IsReactionQuery(self, query):
        """Returns True if this query is likely to be a reaction query.
        
        Args:
            query: the query string.
        """
        m = self.REACTION_MATCHER.match(query.strip())
        return m is not None
    
    def ParseReactionQuery(self, query):
        """Parse the query as a reaction.
        
        Args:
            query: the query string.
        
        Returns:
            An initialized ParsedReaction object, or None if parsing failed.
        """
        try:
            results = self._rparser.parseString(query)
            reactants, products = results
            return ParsedReactionQuery(reactants, products)
        except Exception, e:
            logging.error(e)
        
        return None
        