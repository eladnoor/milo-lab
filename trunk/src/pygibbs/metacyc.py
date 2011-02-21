import os
import logging
import urllib
import re
from toolbox import util, database
import kegg
import json
from collections import deque

def parse_metacyc_file(filename):
    metacyc_file = open(filename, 'r')
    curr_field = ""
    field_map = {}
    line = '#'; 
    line_counter = 0
    entry2fields_map = {}
    while (line):
        #line = str(metacyc_file.readline().rstrip())
        line = metacyc_file.readline().rstrip()
        line_counter += 1
        if (line.startswith(('#', '/')) and not line.startswith('//') ):
            continue

        field = line.split(' - ')[0]
        value = "".join(line.split(' - ')[1:])

        if (field == "//"):
            entry = field_map["UNIQUE-ID"]
            entry2fields_map[entry] = field_map
            field_map = {}
        else:
            if (field != ""):
                curr_field = field
            if (curr_field in field_map):
                field_map[curr_field] = field_map[curr_field] + "\t" + value
            else:
                field_map[curr_field] = value

    
    metacyc_file.close()
    return entry2fields_map

def parse_rxns_metacyc_file(filename):
    metacyc_file = open(filename, 'r')
    prev_field = ""
    prev_value = ""
    field_map = {}
    sparse_map = {}
    line = '#'; 
    line_counter = 0
    entry2fields_map = {}
    non_specific = False
    
    while (line):
        line = str(metacyc_file.readline().rstrip())
        line_counter += 1
        if (line.startswith(('#', '/')) and not line.startswith('//') ):
            continue

        field = line.split(' - ')[0]
        value = "".join(line.split(' - ')[1:])

        if (field == "//"):
            entry = field_map["UNIQUE-ID"]
            field_map['SPARSE'] = dict(filter(lambda k:k[1] != 0, [(k,v) for k,v in sparse_map.iteritems()]))
            if (not non_specific and len(field_map['SPARSE']) > 0):
                entry2fields_map[entry] = field_map
            non_specific = False
            field_map = {}
            sparse_map = {};
        else:
            if (field == 'LEFT' or field == 'RIGHT'):
                prev_field = field
                value = re.sub('\|', '', value)
                prev_value = value
                sparse_map[value] = (sparse_map[value] if (value in sparse_map) else 0) + (1 if (field == 'RIGHT') else -1)
            elif (field.startswith('^COEFFICIENT')):
                try:
                    sparse_map[prev_value] += (int(value) - 1) * (1 if (prev_field == 'RIGHT') else -1) # Value - 1 since 1/-1 was already added in the previous line
                except ValueError:
                    non_specific = True
                    
            elif (field != ""):
                if (field in field_map):
                    field_map[field] = field_map[field] + "\t" + value
                else:
                    field_map[field] = value

    
    metacyc_file.close()
    return entry2fields_map

def parse_metacyc_reaction_formula(formula, uid2compound_map):
    """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
        return the set of substrates, products and the direction of the reaction
    """
    tokens = re.findall("([^=^<]+) (<*=>*) ([^=^>]+)", formula)
    if (len(tokens) != 1):
        raise MetaCycParseException("Cannot parse this formula: " + formula)
    
    (left, direction, right) = tokens[0] # the direction: <=, => or <=>
    
    sparse_reaction = {}
    for (uid, count) in parse_metacyc_reaction_formula_side(left).iteritems():
        comp_uid = re.sub('\|', '', uid)
        if (comp_uid in uid2compound_map):
            sparse_reaction[comp_uid] = sparse_reaction.get(comp_uid, 0) - count
        else:
            raise MetaCycParseException("Compound " + uid + " not in compounds list") 

    for (uid, count) in parse_metacyc_reaction_formula_side(right).iteritems():
        comp_uid = re.sub('\|', '', uid)

        if (comp_uid in uid2compound_map):
            sparse_reaction[comp_uid] = sparse_reaction.get(comp_uid, 0) + count 
        else:
            raise MetaCycParseException("Compound " + uid + " not in compounds list")
 

    return (sparse_reaction, direction)

def parse_metacyc_reaction_formula_side(s):
    """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        return the set of CIDs, ignore stoichiometry
    """
    if (s.strip() == "null"):
        return {}
    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if (len(tokens) == 1):
            amount = 1
            key = member
        else:
            try:
                amount = float(tokens[0])
            except ValueError:
                raise MetaCycParseException("Non-specific reaction: " + s)
            key = tokens[1]
            
        try:
            compound_bag[key] = compound_bag.get(key, 0) + amount
        except ValueError:
            raise MetaCycParseException("Non-specific reaction: " + s)
    
    return compound_bag

def unparse_metacyc_reaction_formula(sparse, direction='=>'):
    s_left = []
    s_right = []
    for (uid, count) in sparse.iteritems():
                
        if (count > 0):
            if (count == 1):
                s_right.append(uid)
            else:
                s_right.append('%d %s' % (count, uid))
        elif (count < 0):
            if (count == -1):
                s_left.append(uid)
            else:
                s_left.append('%d %s' % (-count, uid))
    return ' + '.join(s_left) + ' ' + direction + ' ' + ' + '.join(s_right)

class MetaCycParseException(Exception):
    pass

class MetaCycNonCompoundException(Exception):
    pass

class MetaCyc(object):
    kegg_cids = {} # class static variable
    
    def __init__(self, org='ecoli', db=None):
        self.db = db
        self.org = org
        self.base_dir = '../res/MetaCyc/' + org
        util._mkdir(self.base_dir)
 
        self.TAR_URL = 'http://brg.ai.sri.com/ecocyc/dist/flatfiles-52983746/' + org + '.tar.gz'
        self.TAR_FILE = self.base_dir + '/' + org + '.tar.gz'
       
        self.COMPOUND_FILE = self.base_dir + '/14.6/data/compounds.dat'
        self.REACTION_FILE = self.base_dir + '/14.6/data/reactions.dat'
        self.PATHWAY_FILE = self.base_dir + '/14.6/data/pathways.dat'
        
        if not self.db:
            self.FromFiles()
        elif not self.db.DoesTableExist('metacyc_' + org + '_compound'):
            self.FromFiles()
            self.ToDatabase()
        else:
            self.FromDatabase()
    
    def FromFiles(self):
        self.name2cid_map = {}
        self.uid2compound_map = {}
        self.uid2reaction_map = {}
        self.uid2pathway_map = {}

        logging.info("Retrieving COMPOUNDS file and parsing it")
        if (not os.path.exists(self.COMPOUND_FILE)):
            urllib.urlretrieve(self.TAR_URL, self.TAR_FILE)
            os.chdir(self.base_dir)
            os.system('tar xvfz ' + self.org + '.tar.gz')   

        entry2fields_map = parse_metacyc_file(self.COMPOUND_FILE)
        
        for uid in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[uid]
            
            comp = Compound(uid)
            
            if ("COMMON-NAME" in field_map):
                comp.name = re.sub('<.+?>', '', field_map["COMMON-NAME"].strip())
            if ("SYNONYMS" in field_map):
                all_names = field_map["SYNONYMS"].split('\t')
                for name in all_names:
                    name = re.sub('<.+?>', '', name.strip())
                    self.name2cid_map[name] = uid
                comp.all_names = all_names
            if ("MOLECULAR-WEIGHT" in field_map):
                comp.mass = float(field_map["MOLECULAR-WEIGHT"])
            if ("CHEMICAL-FORMULA" in field_map):    
                comp.formula = field_map["CHEMICAL-FORMULA"]
            if ("INCHI" in field_map):    
                comp.inchi = field_map["INCHI"]
            if ("SMILES" in field_map):    
                comp.smiles = field_map["SMILES"]
            if ("DBLINKS" in field_map):
                for sid in re.findall("PUBCHEM \"(\d+)\"", field_map["DBLINKS"]):
                    comp.pubchem_id = int(sid)
                for cas in re.findall("CAS \"([\d\-]+)\"", field_map["DBLINKS"]):
                    comp.cas = cas
            if ("REGULATES" in field_map):
                comp.regulates = field_map["REGULATES"].split('\t')
            if ("TYPES" in field_map):
                comp.types = field_map["TYPES"].split('\t')
            
            self.uid2compound_map[uid] = comp
        
        logging.info("Retrieving REACTIONS file and parsing it")
        if (not os.path.exists(self.REACTION_FILE)):
            urllib.urlretrieve(self.TAR_URL, self.TAR_FILE)
            os.chdir(self.base_dir)
            os.system('tar xvfz ' + self.org + '.tar.gz')   

        entry2fields_map = parse_rxns_metacyc_file(self.REACTION_FILE)
        
        for uid in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[uid]
            
            direction = '<=>'
            if ("REACTION-DIRECTION" in field_map):
                if (re.search('LEFT-TO-RIGHT', field_map['REACTION-DIRECTION'])):
                    direction = '=>'
                elif (re.search('RIGHT-TO-LEFT', field_map['REACTION-DIRECTION'])):
                    direction = '<='

            rxn = Reaction(uid, sparse_reaction=field_map['SPARSE'], direction=direction)

            if ("COMMON-NAME" in field_map):
                rxn.name = field_map["COMMON-NAME"].strip()
            if ("TYPES" in field_map):
                rxn.types = field_map["TYPES"].split('\t')
            if ("EC-NUMBER" in field_map):    
                rxn.ec_number = field_map["EC-NUMBER"]
            
            self.uid2reaction_map[uid] = rxn

        logging.info("Retrieving PATHWAYS file and parsing it")
        if (not os.path.exists(self.PATHWAY_FILE)):
            urllib.urlretrieve(self.TAR_URL, self.TAR_FILE)
            os.chdir(self.base_dir)
            os.system('tar xvfz ' + self.org + '.tar.gz')   

        entry2fields_map = parse_metacyc_file(self.PATHWAY_FILE)
        
        for uid in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[uid]
            rxn_direction_map = {}
            if ('Super-Pathways' in field_map['TYPES']):
                continue
            
            pw = Pathway(uid)

            if ("COMMON-NAME" in field_map):
                pw.name = field_map["COMMON-NAME"].strip()
            if ("TYPES" in field_map):
                pw.types = field_map["TYPES"].split('\t')
            if ("PREDECESSORS" in field_map):    
                pw.preds = field_map["PREDECESSORS"].split('\t')
                try:
                    pw.UpdateRxnsDict()
                except MetaCycParseException, e:
                    print e
                    continue
            if ("REACTION-LAYOUT" in field_map):
                items = field_map["REACTION-LAYOUT"].split('\t')
                for item in items:
                    tokens = re.findall("\(([^ ]+).+DIRECTION :([^\)]+)", item)
                    if (len(tokens) != 1):
                        raise MetaCycParseException("Pathway %s: Cannot parse reaction layout: %s " % (pw.name, item))
                    (rxn, direction) = tokens[0]
                    if (rxn in rxn_direction_map):
                        raise MetaCycParseException("Pathway %s: Duplicate reaction layout: %s " % (pw.name, item))
                    
                    if (direction == 'L2R'):
                        rxn_direction_map[rxn] = 1
                    elif (direction == 'R2L'):
                        rxn_direction_map[rxn] = -1
                    else:
                        raise MetaCycParseException("Pathway %s: Unknown direction in reaction layout: %s " % (pw.name, item))
                pw.rxn_dirs = rxn_direction_map
                
            self.uid2pathway_map[uid] = pw
                
        
    def ToDatabase(self):
        self.db.CreateTable('metacyc_' + self.org + '_compound', 'uid TEXT, name TEXT, all_names TEXT, '
           'mass REAL, formula TEXT, inchi TEXT, pubchem_id INT, cas TEXT, '
           'regulates TEXT, types TEXT, smiles TEXT')
        self.db
        for uid, comp in self.uid2compound_map.iteritems():
            self.db.Insert('metacyc_' + self.org + '_compound', [uid, comp.name, ';'.join(comp.all_names),
                comp.mass, comp.formula, comp.inchi, comp.pubchem_id, comp.cas,
                ';'.join(comp.regulates), ';'.join(comp.types), comp.smiles])

        self.db.CreateTable('metacyc_' + self.org + '_reaction', 'uid TEXT, name TEXT,  '
                       'ec_number TEXT, types TEXT, equation TEXT')
        for uid, reaction in self.uid2reaction_map.iteritems():
            self.db.Insert('metacyc_' + self.org + '_reaction', [uid, reaction.name, reaction.ec_number, 
                                                    ';'.join(reaction.types), reaction.equation])
        
        self.db.CreateTable('metacyc_' + self.org + '_pathway', 'uid TEXT, name TEXT,  '
                       'types TEXT, preds TEXT, rxn_dirs TEXT')
        for uid, pathway in self.uid2pathway_map.iteritems():
            self.db.Insert('metacyc_' + self.org + '_pathway', [uid, pathway.name, ';'.join(pathway.types), ';'.join(pathway.preds), json.dumps(pathway.rxn_dirs)])
            
        self.db.Commit()


    def FromDatabase(self):
        logging.info('Reading ' + self.org + ' MetaCyc data from the database')
        
        self.name2cid_map = {}
        self.uid2compound_map = {}
        self.uid2reaction_map = {}
        self.uid2pathway_map = {}

        for row in self.db.DictReader('metacyc_' + self.org + '_compound'):
            comp = Compound(row['uid'], row['name'], row['all_names'].split(';'), 
                            row['mass'], row['formula'], row['inchi'], row['pubchem_id'],
                            row['cas'], row['regulates'].split(';'), row['types'].split(';'), row['smiles'])
            self.uid2compound_map[row['uid']] = comp
        
        for row in self.db.DictReader('metacyc_' + self.org + '_reaction'):
            try:
                (sparse, direction) = parse_metacyc_reaction_formula(row['equation'], self.uid2compound_map)
                reaction = Reaction(uid=row['uid'], name=row['name'], sparse_reaction=sparse, 
                                 direction=direction, types=row['types'].split(';'), ec_number=row['ec_number'])
                self.uid2reaction_map[row['uid']] = reaction
            except MetaCycParseException, e:
                print "Cannot parse equation: (%s) for reaction: %s: %s" % (row['equation'], row['uid'], e)

        for row in self.db.DictReader('metacyc_' + self.org + '_pathway'):
            try:
                pathway = Pathway(uid=row['uid'], name=row['name'], types=row['types'].split(';'), 
                                  preds=row['preds'].split(';'), rxn_dirs=json.loads(row['rxn_dirs']))
                self.uid2pathway_map[row['uid']] = pathway
            except MetaCycParseException, e:
                print e
    
    def rxn_uid2sparse_reaction(self, uid):
        if (uid in self.uid2reaction_map):
            return self.uid2reaction_map[uid].sparse
        else:
            return None
    
    def sparse2kegg_cids(self, sparse, kegg_inst):
        if (len(self.kegg_cids) == 0):
            for name,cid in kegg_inst.name2cid_map.iteritems():
                self.kegg_cids[name.lower()] = cid

        kegg_dict = {}
     
        for (uid, value) in sparse.iteritems():
            name_cand_cid = None
            inchi_cand = None
            if (not uid in self.uid2compound_map):
                raise MetaCycNonCompoundException('Compound %s not found' % uid)
            comp = self.uid2compound_map[uid]
            if (comp.name.lower() in self.kegg_cids):
                name_cand_cid = self.kegg_cids[comp.name.lower()] 
            else:
                for syn in comp.all_names:
                    if (syn.lower() in self.kegg_cids):
                        name_cand_cid = self.kegg_cids[syn.lower()]
                        break
            inchi = comp.inchi
            if (inchi):
                inchi_cand = kegg_inst.inchi2cid(inchi) 
            
            if (name_cand_cid and inchi_cand and name_cand_cid != inchi_cand):
                logging.info('Ambigious mapping of MetaCyc compound %s to Kegg, matches cid %s by name and cid %s by inchi, selecting inchi...' % (uid, name_cand_cid, inchi_cand))
                kegg_dict[inchi_cand] = value 
            elif (not name_cand_cid and not inchi_cand):
                raise kegg.KeggNonCompoundException('No mapping of MetaCyc compound %s to Kegg' % uid)
            else:
                kegg_dict[name_cand_cid or inchi_cand] = value
                
        return kegg_dict

class Compound(object):
        
    def __init__(self, uid=None, name=None, all_names=None, mass=None,
                 formula=None, inchi=None, pubchem_id=None, cas=None, 
                 regulates=None, types=None, smiles=None):

            self.uid = uid;                     # UNIQUE-ID
            self.name = name                    # COMMON-NAME
            if (self.name):
                self.name = re.sub('<.+?>', '', name) # Removing HTML tags
            self.all_names = []                 # SYNONYMS
            if (all_names and len(all_names) > 0):
                for s in all_names:
                    self.all_names.append(re.sub('<.+?>', '', s))
            self.mass = mass                    # MOLECULAR-WEIGHT
            self.formula = formula              # CHEMICAL-FORMULA
            self.inchi = inchi if inchi != None else ""                  # INCHI
            self.pubchem_id = None              # Parsed from DBLINKS
            self.cas = ""                       # Parsed from DBLINKS
            self.regulates = regulates if regulates != None else []                 # REGULATES
            self.types = types if types != None else []                     # TYPES
            self.smiles = smiles if smiles != None else ""                 # SMILES
            if (smiles and not inchi):
                self.inchi = kegg.smiles2inchi(smiles)
            
    def Print(self):
        print 'uid=',self.uid,'\nname=',self.name,'\nall_names=',self.all_names,'\nmass=',self.mass,'\nformula=',self.formula,'\ninchi=',self.inchi,'\npubchem_id=',self.pubchem_id,'ncas=',self.cas,'\nregs=',self.regulates,'\ntypes=',self.types,'\nsmiles=',self.smiles,'\n' 
        
class Reaction(object):
    
    def __init__(self, uid=None, name=None, sparse_reaction=None, direction=None, types=None, ec_number=None):
        self.uid = uid
        self.name = name
        self.sparse = sparse_reaction
        self.direction = direction
        self.types = types;
        self.ec_number = ec_number
        self.equation = unparse_metacyc_reaction_formula(sparse_reaction, direction)   
    
class Pathway(object):
    
    def __init__(self, uid=None, name=None, types=None, preds=None, rxn_dirs=None):
        self.uid = uid
        self.name = name
        self.types = types if types != None else [] 
        self.preds = preds if preds != None else []
        self.rxns_dict = self.Preds2RxnsOrder() if preds != None else {}
        self.rxn_dirs = rxn_dirs if rxn_dirs != None else {}
        
    def GetRxnsOrder(self):
        return self.rxns_dict
    
    def GetRxnsDirs(self):
        return self.rxn_dirs
    
    def UpdateRxnsDict(self):
        self.rxns_dict = self.Preds2RxnsOrder() if self.preds != None else {}
        
    def Preds2RxnsOrder(self):
        rxn2position = {}
        pw_parent = ''
        if (len(self.preds) > 0):
            pairs_map = {}
            parents = []
            sons = []
            
            # Loading pairs to map
            for pair in self.preds:
                vals =  re.sub('[\(\)]', '', pair).split(' ')
                if (len(vals) == 1): # Single reaction - no parents. Update children positions
                    if (pw_parent):
                        raise  MetaCycParseException("More then single starting reaction in pathway " + self.name)
                    pw_parent = vals[0]
                elif (len(vals) == 2):
                    son,parent = vals[0],vals[1]
                    parents.append(parent)
                    sons.append(son)
                    childs = pairs_map.get(parent, '')
                    if (len(childs) == 0):
                        pairs_map[parent] = son
                    else:
                        pairs_map[parent] = childs + ';' + son
                        
            if (len(pw_parent) == 0):
                parents_cands = set(parents) - set(sons)
                if (len(parents_cands) != 1):
                    raise MetaCycParseException("Expecting single parent, found %d in pathway: %s" % (len(parents_cands), self.name))
                else:
                    pw_parent = parents_cands.pop()
            
            rxn2position[pw_parent] = 1
            curr_parent = pw_parent
            parents_q = deque([curr_parent])
            
            while (len(parents_q) > 0):
                curr_parent = parents_q.popleft()
                if (curr_parent in pairs_map):
                    for child in pairs_map[curr_parent].split(';'):
                        if (child in rxn2position):
                            logging.info("Circular pathway: " + self.name)
                        else:
                            rxn2position[child] = rxn2position[curr_parent] + 1
                            parents_q.append(child)
        return rxn2position
    

                
                
##########################
##  ##  ##  ##  ##  ##  ##
##  ##  ##  ##  ##  ##  ##
##########################     
def test():
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    #db = database.SqliteDatabase('../MetaCyc/yanivtest.sqlite')
    kg = kegg.Kegg.getInstance()
    ecocyc = MetaCyc('ecoli', db)
    comps = ecocyc.uid2compound_map
    rxns = ecocyc.uid2reaction_map
    pws = ecocyc.uid2pathway_map
    print '#comps= ',comps.__len__(), '\n#reactions= ', rxns.__len__(),'\n#pathways= ',pws.__len__()
    
    
    
        
if __name__ == "__main__":
    test()
 
    