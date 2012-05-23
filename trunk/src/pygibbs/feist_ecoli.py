import re, csv
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_errors import KeggReactionNotBalancedException
import logging
from pygibbs.kegg import Kegg
import pybel
import openbabel
from pygibbs.thermodynamic_estimators import LoadAllEstimators

class IrrevParseException(Exception):
    pass

class Feist():
    
    @staticmethod
    def NormalizeInChI(inchi):
        
        # TODO: 
        # There are still some mapping problems that are caused by discarding all
        # these layers:
        # 1) Fe2+ and Fe3+ are mapped both to C00023, although they have a different number of electrons
        # 2) The InChI for heme does not include the metal ion or something
        
        if not inchi:
            return None
        for layer in ['p', 'm', 's', 't', 'b', 'q']:
            inchi = re.sub('/' + layer + '([0-9\+\-\;\?\,]+)', '', inchi)
        #inchi = re.sub('/h([0-9\+\-\;\?\,\(\)H]+)', '', inchi)
        return inchi
    
    @staticmethod
    def ReadBiGGCompounds():
        sdf_input_filename = '../data/metabolic_models/iAF1260.sdf'
        sdfile = pybel.readfile("sdf", sdf_input_filename)
        biggID2inchi = {}
        
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("inchi")
        for m in sdfile:
            inchi = obConversion.WriteString(m.OBMol).strip()
            inchi = Feist.NormalizeInChI(inchi)
            if not inchi:
                inchi = None
            biggID = m.title.strip()
            biggID2inchi[biggID] = inchi
        return biggID2inchi
    
    @staticmethod
    def ReadKeggCompounds():
        kegg = Kegg.getInstance()
        inchi2KeggID = {}
        inchi2KeggID[None] = 0
        for cid in sorted(kegg.get_all_cids()):
            inchi = kegg.cid2inchi(cid)
            inchi = Feist.NormalizeInChI(inchi)
            if inchi not in inchi2KeggID: # since CIDs are sorted, this will always keep the lowest CID with this InChI
                inchi2KeggID[inchi] = cid
        return inchi2KeggID

    @staticmethod
    def Bigg2KEGG():
        bigg2inchi = Feist.ReadBiGGCompounds()
        inchi2kegg = Feist.ReadKeggCompounds()

        bigg2kegg = {}
        for biggID in bigg2inchi.keys():
            bigg2kegg[biggID] = inchi2kegg.get(bigg2inchi[biggID], 0)
        
        # extras
#        bigg2kegg['h'] = 80
#        bigg2kegg['adocbip'] = 6509
#        bigg2kegg['agdpcbi'] = 6510
#        bigg2kegg['adocbi'] = 6508
#        bigg2kegg['adocbl'] = 0
#        bigg2kegg['s'] = 87
#        bigg2kegg['cbi'] = 5774
#        bigg2kegg['cbl1'] = 0
#        bigg2kegg['pppg9'] = 1079
#        bigg2kegg['ppp9'] = 2191
#        bigg2kegg['pheme'] = 32
#        bigg2kegg['dsbard'] = 0
#        bigg2kegg['dsbaox'] = 0
#        bigg2kegg['dsbcox'] = 0
#        bigg2kegg['dsbcrd'] = 0
#        bigg2kegg['dsbdox'] = 0
#        bigg2kegg['dsbdrd'] = 0
#        bigg2kegg['dsbgox'] = 0
#        bigg2kegg['dsbgrd'] = 0
#        bigg2kegg['fldox'] = 2869
#        bigg2kegg['fldrd'] = 2745
#        bigg2kegg['hemeO' ] = 0
#        bigg2kegg['no'] = 533
        

        return bigg2kegg

    @staticmethod
    def _Bigg2KEGG():
        compound_file = '../data/metabolic_models/iAF1260_compounds.csv'
        bigg2kegg = {}
        for row in csv.DictReader(open(compound_file, 'r')):
            if row['KeggID']:
                bigg2kegg[row['abbreviation']] = int(row['KeggID'][1:])
            else:
                bigg2kegg[row['abbreviation']] = 0
        return bigg2kegg

    def __init__(self):
        self.reactions = []
        self.bigg2kegg = Feist.Bigg2KEGG()
        for biggID, keggID in Feist._Bigg2KEGG().iteritems():
            if biggID not in self.bigg2kegg or self.bigg2kegg[biggID] == 0:
                self.bigg2kegg[biggID] = keggID
    
    @staticmethod
    def FromFiles():
        feist = Feist()
        reaction_file = '../data/metabolic_models/iAF1260_reactions.csv'
        
        # Read compounds ids 2 Kegg cid's mapping file into a dict
    
        counters = {'kegg_error':0, 'translocation':0, 'exchange':0, 'sink':0,
                    'unbalanced':0, 'okay':0}
        for row in csv.DictReader(open(reaction_file, 'r')):
            #if 'Transport' in row['subSystem']:
            #    counters['translocation'] += 1
            #    continue
            if row['abbreviation'][0:3] == 'EX_':
                counters['exchange'] += 1
                continue
            if row['abbreviation'][0:3] == 'DM_':
                counters['sink'] += 1
                continue            
            
            sparse = Feist.parse_equation(row['equation'])
            kegg_sparse = {}
            for biggID, coeff in sparse.iteritems():
                keggID = feist.bigg2kegg[biggID]
                kegg_sparse[keggID] = kegg_sparse.get(keggID, 0) + coeff

            if 0 in kegg_sparse:
                logging.debug('Some compounds are missing KEGG IDs in %s: %s' %
                              (row['abbreviation'], row['equation']))
                counters['kegg_error'] += 1
                continue                
                
            for keggID in [k for k, v in kegg_sparse.iteritems() if v == 0]:
                del kegg_sparse[keggID]
            
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
                reaction.Balance(balance_water=True, exception_if_unknown=False)
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
    def parse_equation(s_input):
        equation = re.sub('\[[pce]\]', '', s_input)
        equation = re.sub(' : ', '', equation)
        left, right = re.split('<==>|-->', equation)
        
        sparse = {}
        if left:
            Feist.parse_equation_side(left, -1, sparse)
        if right:
            Feist.parse_equation_side(right, 1, sparse)
        return sparse

    @staticmethod
    def parse_equation_side(s_input, coeff, sparse):
        for comp in s_input.split('+'):
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
    
    estimators = LoadAllEstimators()
    thermo = estimators['UGC']
    thermo.SetConditions(pH=7, I=0.15, T=298.15, pMg=14)
    
    csv_out = csv.writer(open('../res/iAF1260_thermo.csv', 'w'))
    csv_out.writerow(['reaction', 'dG0_prime'])
    dG0s = thermo.GetTransfromedKeggReactionEnergies(feist.reactions)
    for i, reaction in enumerate(feist.reactions):
        csv_out.writerow([reaction.name, dG0s[0, i]])
        
        
        
    