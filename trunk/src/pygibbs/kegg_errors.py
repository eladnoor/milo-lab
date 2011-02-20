#!/usr/bin/python

class KeggParseException(Exception):
    pass
        
class KeggNonCompoundException(Exception):
    pass

class KeggReactionNotBalancedException(Exception):
    pass
    
class KeggMissingModuleException(Exception):
    pass