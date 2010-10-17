class Compound(object):
    """A representation of a single compound."""    
    def __init__(self, name, kegg_id, formula):
        """Initialize a compound."""
        self.kegg_id = kegg_id
        self.formula = formula
        self.name = name
        
    def HasData(self):
        return self.kegg_id or self.formula or self.name
    
    def IsComplete(self):
        return self.kegg_id and self.formula and self.name
        
    def __str__(self):
        return self.name

def ReadCompoundsFromCsvLine(line):
    """Reads a list of compounds from my CSV file.
    
    Arguments:
        line: the line of CSV to read from.
        
    Returns:
        A dictionary mapping compound names to compounds.
    """
    fields = line.strip().split(', ')
    kegg_id = fields[0]
    formula = fields[1]
    
    compounds = {}
    for name in fields[2:]:
        c = Compound(name, kegg_id, formula)
        if c.IsComplete():
            compounds[name] = c
    return compounds


def ReadCompoundsFromCsvFile(filename):
    compounds = {}
    f = open(filename, 'r')
    for line in f:
        compounds.update(ReadCompoundsFromCsvLine(line))
    
    f.close()
    return compounds
    
