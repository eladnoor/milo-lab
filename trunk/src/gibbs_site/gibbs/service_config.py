from gibbs import models
from matching import approximate_matcher
from matching import reaction_parser


def Singleton(cls):
    """Decorator that makes a class a singleton.
    
    TODO(flamholz): move somewhere more central.
    """
    instance_container = []
    def getinstance():
        if not len(instance_container):
            instance_container.append(cls())
        return instance_container[0]
    return getinstance


@Singleton
class ServiceConfig(object):
    """A singleton class that contains global service configuration and state.
    
    All state/config is considered unmodifiable.
    """
    
    def __init__(self):
        self._compound_library = models.Compound.AsLibrary()
        self._compound_matcher = approximate_matcher.BackfillingRegexApproxMatcher(
            self._compound_library, max_results=10, min_score=0.1)
        self._reaction_parser = reaction_parser.ReactionParser(self._compound_matcher)
    
    compound_matcher = property(lambda self: self._compound_matcher)
    reaction_parser = property(lambda self: self._reaction_parser)


def Get():
    """Convenience method to get the single ServiceConfig instance."""
    return ServiceConfig()  # Singleton decorator ensures there's only one