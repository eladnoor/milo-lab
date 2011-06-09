class MissingCompoundFormationEnergy(Exception):
    def __init__(self, value, cid=0):
        self.value = value
        self.cid = cid
    def __str__(self):
        return repr(self.value)
    
class MissingReactionEnergy(Exception):
    def __init__(self, value, sparse):
        self.value = value
        self.sparse = sparse
    def __str__(self):
        return repr(self.value)