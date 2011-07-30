import re, csv
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggNonCompoundException,\
    KeggReactionNotBalancedException
import logging

class IrrevParseException(Exception):
    pass

class Feist():
    
    def __init__(self):
        self.reactions = []
        
    @staticmethod
    def FromFiles():
        feist = Feist()
        comp2cid_file = '../data/metabolic_models/ecoli2kegg_cid.txt'
        rxns_file = '../data/metabolic_models/reactionList_iAF1260.txt'
        
        # Read compounds ids 2 Kegg cid's mapping file into a dict
        comp2cid_map = {}
        for row in csv.reader(open(comp2cid_file, 'r'), delimiter='\t'):
            if len(row) != 2:
                raise Exception("syntax error in %s: %s" % (comp2cid_file, str(row)))
            id, cid = row
            comp2cid_map[id] = cid
    
        rxns_input = open(rxns_file, 'r')
        line = rxns_input.readline().strip()
        
        while line:
            line = rxns_input.readline().strip() # Skipping the header line
            if line == '':
                continue
            
            fields = line.split('\t')
            name = fields[1]
            equation = fields[3]
            reversible = fields[7]
            
            # All fields: abbrev, name, syn, equation, subsystem, compartment, ecnumber, reversible, translocation, internal_id, confidence_score, notes 
    
            equation = re.sub('\[[pce]\]', '', equation)
            equation = re.sub(' : ', '', equation)
            left, right = re.split('<==>|-->', equation)
            
            try:
                sparse = {}
                Feist.parse_palsson_side_eq(left, -1, sparse, comp2cid_map)
                Feist.parse_palsson_side_eq(right, 1, sparse, comp2cid_map)
            except IrrevParseException:
                logging.debug('Error parsing reaction %s: %s' % (name, equation))
                continue
            except KeggNonCompoundException:
                logging.debug('No Kegg cid for reaction %s: %s' % (name, equation))
                continue
            
            for cid in sparse.keys():
                if sparse[cid] == 0:
                    del sparse[cid]
            
            if sparse == {}:
                logging.debug('Reaction is empty %s: %s' % (name, equation))
                continue
            
            if reversible == 'Reversible':
                direction = '<=>'
            elif reversible == 'Irreversible':
                direction = '=>'
            else:
                raise ValueError('unknown irreversibility tag: ' + reversible)
            
            reaction = Reaction(name, sparse, direction=direction)
            try:
                reaction.Balance(balance_water=True, exception_if_unknown=True)
            except KeggReactionNotBalancedException:
                pass
            
            feist.reactions.append(reaction)
        return feist

    @staticmethod
    def parse_palsson_side_eq(str, coef, sparse, comp2cid_map):
        for comp in str.split('+'):
            comp = comp.strip()
            elems = re.findall('(\(.+\))* *(\w+)', comp)
            if (len(elems) != 1):
                raise IrrevParseException
            
            s,comp = elems[0]
            if len(s) == 0:
                s = 1
            else:
                s = float(s[1:-1]) # Removing the parenthesis        
            if (comp in comp2cid_map):
                cid = comp2cid_map[comp]
                if int(cid[1:]) in sparse:
                    sparse[int(cid[1:])] += coef * s
                else:
                    sparse[int(cid[1:])] = coef * s
            else:
                raise KeggNonCompoundException

if __name__ == "__main__":
    #logging.getLogger('').setLevel(logging.DEBUG)
    feist = Feist.FromFiles()
    for reaction in feist.reactions:
        print str(reaction), reaction.direction, str(reaction.sparse)