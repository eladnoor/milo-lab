import re, csv
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggReactionNotBalancedException
import logging

class IrrevParseException(Exception):
    pass

class Feist():
    
    kegg_misses = {}
    
    def __init__(self):
        self.reactions = []
        
    @staticmethod
    def FromFiles():
        feist = Feist()
        compound_file = '../data/metabolic_models/iAF1260_compounds.csv'
        reaction_file = '../data/metabolic_models/iAF1260_reactions.csv'
        
        # Read compounds ids 2 Kegg cid's mapping file into a dict
        comp2cid_map = {}
        for row in csv.DictReader(open(compound_file, 'r')):
            if row['KeggID']:
                comp2cid_map[row['abbreviation']] = int(row['KeggID'][1:])
            else:
                comp2cid_map[row['abbreviation']] = 0
    
        counters = {'kegg_error':0, 'translocation':0, 'exchange':0, 'sink':0,
                    'unbalanced':0, 'okay':0}
        for row in csv.DictReader(open(reaction_file, 'r')):
            if 'Transport' in row['subSystem']:
                counters['translocation'] += 1
                continue
            if row['abbreviation'][0:3] == 'EX_':
                counters['exchange'] += 1
                continue
            if row['abbreviation'][0:3] == 'DM_':
                counters['sink'] += 1
                continue            
            
            sparse = Feist.parse_equation(row['equation'])
            kegg_sparse = dict([(comp2cid_map[comp], coef) for 
                                (comp, coef) in sparse.iteritems()])
            if 0 in kegg_sparse:
                logging.debug('Some compounds are missing KEGG IDs in %s: %s' %
                              (row['abbreviation'], row['equation']))
                counters['kegg_error'] += 1
                continue                
            
            directionality = row["directionality without uncertainty (pH 7.2)"]
            #directionaliry = row['reconstruction directionality']
            if directionality == 'reversible':
                direction = '<=>'
            elif directionality == 'forward only':
                direction = '=>'
            elif directionality == 'reverse only':
                direction = '<='
            else:
                raise ValueError('unknown directionality tag: ' + directionality)
            
            reaction = Reaction(row['abbreviation'], kegg_sparse, direction=direction)
            try:
                reaction.Balance(balance_water=True, exception_if_unknown=True)
                counters['okay'] += 1
            except KeggReactionNotBalancedException as e:
                logging.debug(str(e) + ' - ' + str(reaction))
                counters['unbalanced'] += 1
                continue
            
            feist.reactions.append(reaction)
        
        logging.debug(" ; ".join(["%s : %d" % (key, val)
                                for (key, val) in counters.iteritems()]))
        return feist

    @staticmethod
    def parse_equation(str):
        equation = re.sub('\[[pce]\]', '', str)
        equation = re.sub(' : ', '', equation)
        left, right = re.split('<==>|-->', equation)
        
        sparse = {}
        if left:
            Feist.parse_equation_side(left, -1, sparse)
        if right:
            Feist.parse_equation_side(right, 1, sparse)
        return sparse

    @staticmethod
    def parse_equation_side(str, coeff, sparse):
        for comp in str.split('+'):
            comp = comp.strip()
            elems = re.findall('(\(.+\))* *([\w\-\(\)]+)', comp)
            if len(elems) != 1:
                raise IrrevParseException
            
            s, comp = elems[0]
            if len(s) == 0:
                sparse[comp] = sparse.get(comp, 0) + coeff
            else:
                s = float(s[1:-1]) # Removing the parentheses        
                sparse[comp] = sparse.get(comp, 0) + coeff * s

if __name__ == "__main__":
    logging.getLogger('').setLevel(logging.DEBUG)
    feist = Feist.FromFiles()
    #for reaction in feist.reactions:
    #    print str(reaction)