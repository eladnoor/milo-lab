import re
import csv
import numpy as np
from pygibbs.thermodynamic_constants import default_pH
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
import types

BUFFERS_CSV_FNAME = '../data/thermodynamics/pKa_of_buffers.csv'
NIST_ADDITIONS_FILE = '../data/thermodynamics/nist_additions.csv'
missing_set = set()

class MissingCompoundsFromKeggException(Exception):
    def __init__(self, names):
        self.names = names
    def __str__(self):
        return "Could not find this reactant in the list: " + \
            ', '.join(self.names)

def remove_superfluous_chars(s):
    s = s.replace('"', '')
    s = s.strip()
    if s[0] == '-':
        s = s[1:]
    return s.lower()

def parse_single_reactant(s, compound_aliases):
    s = remove_superfluous_chars(s)

    # Water is the only compound allowed to be in liquid phase
    if s[-6:] == 'h2o(l)':
        s = s[:-3]

    # cannot use compounds in gas phase or non-aqueous solvents
    if s[-3:] in ['(g)', '(l)', '(s)'] or s[-5:] == '(sln)':
        raise MissingCompoundsFromKeggException([s])

    if s[-4:] == '(aq)':
        s = s[:-4]
    
    if s == "1/2 o2": # special case that confuses the regular expression
        return "0.5 C00007"
    if s in compound_aliases:
        return compound_aliases[s]

    tmp = re.findall('^(\d+) (.*)', s)
    if tmp:
        count, name = tmp[0]
        if name in compound_aliases:
            if compound_aliases[name] is None:
                return None
            else:
                return count + " " + compound_aliases[name]
    
    raise MissingCompoundsFromKeggException([s])                

def parse_reaction_side(s, compound_aliases):
    res = []
    missing_names = []
    for r in s.split(' + '):
        try:
            reactant = parse_single_reactant(r, compound_aliases)
            if reactant is not None:
                res.append(reactant)
        except MissingCompoundsFromKeggException as e:
            missing_names += e.names
    if missing_names:
        raise MissingCompoundsFromKeggException(missing_names)
    else:
        return " + ".join(res)

def parse_reaction(s, compound_aliases):
    if s.find('=') == -1:
        raise MissingCompoundsFromKeggException([])
    
    subs, prods = s.split('=', 1)
    missing_names = []
    try:
        subs = parse_reaction_side(subs, compound_aliases)
    except MissingCompoundsFromKeggException as e:
        missing_names += e.names

    try:
        prods = parse_reaction_side(prods, compound_aliases)
    except MissingCompoundsFromKeggException as e:
        missing_names += e.names

    if missing_names:
        raise MissingCompoundsFromKeggException(missing_names)
    else:
        return subs + ' = ' + prods
    
def make_sure_it_is_float(s):
    if not s or s == 'not given' or s == '-':
        return None
    
    if s == '298.l5':
        return 298.15

    if s[0] == '~':
        s = s[1:]

    if s[0:2] == '? ':
        s = s[2:]

    if s.find('?') != -1:
        return None

    tmp = re.findall('(\d+\.?\d*) \- (\d+\.?\d*)', s)
    if tmp:
        v1, v2 = tmp[0]
        v = (float(v1) + float(v2)) / 2.0
        #print "%s -> %g" % (s, v)
        return v

    tmp = re.findall('(\d+\.?\d*) 10(\-?\d+)', s)
    if tmp:
        v1, v2 = tmp[0]
        v = float(v1) * 10**float(v2)
        #print "%s -> %g" % (s, v)
        return v

    return float(s)

def choose_buffer(row_dict):
    """ choose the first buffer that has a concentration value """
    
    for key in ['buffer', 'buffer(mol dm-3)', 'Buffer']:
        if row_dict[key] and re.findall('\d+\.\d+', row_dict[key]):
            return row_dict[key].lower()
    return None

def get_buffer_charges(base_charge, pKa_list, conc, pH):
    if not pH:
        pH = default_pH
    
    pKa_list = sorted(pKa_list, reverse=True)
    
    species_proportions = []
    for n in xrange(len(pKa_list) + 1):
        p = np.prod([10**(pKa - pH) for pKa in pKa_list[:n]])
        species_proportions.append(p)
    
    total = sum(species_proportions)
    species_concentration = [(conc * p / total) for p in reversed(species_proportions)]
    # Note that the calculation of the species lists them in order of increasing
    # charges. Since we are given the base-charge (i.e. the highest charge) we
    # need to reverse the order to so that the first value will correspond to
    # the bast_charge 
    
    charge_conc_pairs = []
    for i, c in enumerate(species_concentration):
        charge_conc_pairs.append((base_charge - i, c))
    return charge_conc_pairs

def buffer_match(buffer_name, buffer_dict, pH, missing_buffers):
    if not buffer_name:
        return []
    
    buffer_name = buffer_name.replace('_-', '')
    buffer_name = buffer_name.replace('{', '')
    buffer_name = buffer_name.replace('}', '')

    if buffer_name.find(' or ') != -1:
        for b in re.split(' or ', buffer_name):
            return buffer_match(b, buffer_dict, pH, missing_buffers)

    sub_buffer_names = re.split(' and |\s?\+ |, |\/| and\/or ', buffer_name)
    if len(sub_buffer_names) > 1:
        res = []
        for b in sub_buffer_names:
            res += buffer_match(b, buffer_dict, pH, missing_buffers)
        return res
    
    tmp = re.findall('^(\w+)\s?.?(\d+\.\d+)( m)*(\))?$', buffer_name)
    if tmp:
        buffer_name = "%s (%s mol dm-3)" % tmp[0][0:2]
    buffer_name.replace('kh2po4', 'potassium phosphate')
    
    tmp = re.findall('^([0-9\s\w\-\,\[\]\(\)]+)\((.?\d+\.\d+.?)\smol (dm-3)?(kg-1)?\) *(\+ hcl)?(\+ koh)?(\+ naoh)?$', buffer_name)
    if tmp:
        b = tmp[0][0].strip()
        conc = make_sure_it_is_float(tmp[0][1])
        if not conc:
            return []

        charge_conc_pairs = []
        if b in buffer_dict:
            base_charge, pKa_list = buffer_dict[b]
            charge_conc_pairs += get_buffer_charges(base_charge, pKa_list, conc, pH)
        else:
            missing_buffers[b] = missing_buffers.get(b, 0) + 1
        
        return charge_conc_pairs

    tmp = re.findall('^(\w+) (\w+)\s*\(.?(\d+\.\d+).?\smol dm-3\) *(\+ hcl)?$', buffer_name)
    if tmp and tmp[0][0] in ['potassium', 'sodium']:
        conc = make_sure_it_is_float(tmp[0][2])
        if not conc:
            return []

        charge_conc_pairs = [(1, conc)] # add the +1 ion for the potassium/sodium

        if re.findall('\+ hcl$', buffer_name):
            charge_conc_pairs.append((-1, conc))

        b = "%s (%s mol dm-3)" % (tmp[0][1], tmp[0][2])
        return charge_conc_pairs + buffer_match(b, buffer_dict, pH, missing_buffers)
    
    b = "<" + buffer_name + ">"
    missing_buffers[b] = missing_buffers.get(b, 0) + 1
    return []
        
def load_buffer_dict():
    """ load the pKa of all the buffers """
    buffer_dict = {}
    for row_dict in csv.DictReader(open(BUFFERS_CSV_FNAME, 'r')):
        buffer_names = row_dict['Buffer']
        base_charge = int(row_dict['base charge'])
        pKa_list = []
        for i in range(1, 5):
            if row_dict.get('pKa%d' % i, None):
                pKa_list.append(float(row_dict['pKa%d' % i]))
        for buffer_name in buffer_names.split(';'):
            buffer_dict[buffer_name] = (base_charge, pKa_list)
    return buffer_dict

def load_compound_aliases():
    """
        create a hash that contains # keys (representative names) and values 
        (array of all known aliases for each representative name)
    """
    compound_aliases = {}
    kegg = Kegg.getInstance()
    for cid in sorted(kegg.cid2compound_map.keys()):
        for alias in kegg.cid2compound_map[cid].all_names:
            alias = remove_superfluous_chars(alias)
            compound_aliases.setdefault(alias, "C%05d" % cid)
    
    compound_aliases['carbon dioxide'] = 'C00288'
    compound_aliases['h2o'] = 'C00001'
    compound_aliases['H+'] = None
    return compound_aliases

################################################################################

def WriteDataToDB(db):
    """
        each reaction is composed of 'substrates' and 'products'.
        each compound from both classes is matched to aliases and receive its 
        representative name
    """
    buffer_dict = load_buffer_dict()
    compound_aliases = load_compound_aliases()
    titles2colnum = {}
    for new_row_dict in db.DictReader('nist_fields'):
        titles2colnum[new_row_dict['name']] = int(new_row_dict['col_number'])
    
    salt_titles = { # values are lists of (charge, conc) pairs
        'm(KCl,mol.kg-3)':    [(1, 1), (-1, 1)], # K(+) Cl(-)
        'c(KCl,mol dm-3)':    [(1, 1), (-1, 1)], # K(+) Cl(-)
        'c(MnCl2,mol dm-3)':  [(2, 1), (-1, 2)], # Mn(2+) 2xCl(-)
        'c(MgCl2)':           [(2, 1), (-1, 2)], # Mg(2+) 2xCl(-)
        'c(MnSO4,mol dm-3)':  [(2, 1), (-2, 1)], # Mn(2+) SO4(2-)
        'c c(Mg2+,mol dm-3)': [(2, 1)],          # Mg(2+)
        'c(Mg)tot(mol dm-3)': [(2, 1)],          # Mg(2+)
        'c(MgSO4,mol dm-3)':  [(2, 1), (-2, 1)], # Mg(2+) SO4(2-)
        'm(MgCl2,mol.kg-1)':  [(2, 1), (-1, 2)], # Mg(2+) 2xCl(-)
        'c(CaCl2,mol dm-3)':  [(2, 1), (-1, 2)], # Ca(2+) 2xCl(-)
        'c(ZnCl2,mol dm-3)':  [(2, 1), (-1, 2)], # Zn(2+) 2xCl(-)
        'm(KCl,mol.kg-1)':    [(1, 1), (-1, 1)], # K(+) Cl(-)
        'c(Na+,mol dm-3)':    [(1, 1)],          # Na(+)
        'c(NaCl,mol dm-3)':   [(1, 1), (-1, 1)], # Na(+) Cl(-)
        'c(MgCl2,mol dm-3)':  [(2, 1), (-1, 2)], # Mg(2+) 2xCl(-)
        'c(Mg2+,mol dm-3)':   [(2, 1)],          # Mg(2+)
        'c(MnSO4)':           [(2, 1), (-2, 1)], # Mn(2+) SO4(2-)
        'pMn':                [(2, 1)], # the logscale is taken care of 
        'pMg':                [(2, 1)], # specifically in the code itself
        'c(orthophosphate)': 'phosphate',
        'c(orthophosphate,mol dm-3)': 'phosphate',
        'c(phosphate,mol dm-3)': 'phosphate',
        'c(Tris,mol dm-3)': 'tris'
        }

    Mg_titles = ['c(MgCl2)', 'c c(Mg2+,mol dm-3)', 'c(Mg)tot(mol dm-3)',
        'c(MgSO4,mol dm-3)', 'm(MgCl2,mol.kg-1)', 'c(MgCl2,mol dm-3)',
        'c(Mg2+,mol dm-3)', 'pMg']
    
    title_mapping = {'URL':'url', 'Reference_id':'reference_id', 
        'Method':'method', 'Evaluation':'evaluation', 'EC value':'ec', 
        'Enzyme':'enzyme', 'Reaction':'reaction', 'pH':'pH',
        'Ic(mol dm-3)':'I', 'Ic':'I', 'Im(mol kg-1)':'I', 'Ic(kJ.mol-1)':'I',
        'T(T)':'T', 'T(K)':'T',
        "K'":"K_tag", "Kc'(mol dm-3)":"K_tag", "Kc'":"K_tag", "K '":"K_tag", "Km'":"K_tag",
        "K'(kJ.mol-1)":"K_tag", "Kc":"K", "K":"K"}
    
    text_titles = ['url', 'reference_id', 
        'method', 'evaluation', 'ec', 'enzyme',
        'kegg_reaction', 'reaction']
    real_titles = ['K', 'K_tag', 'T', 'I', 
        'pH', 'pMg']
    new_titles = text_titles + real_titles
    
    db.CreateTable('nist_equilibrium', [t + " TEXT" for t in text_titles] + [t + " REAL" for t in real_titles])
    db.CreateTable('nist_errors', ['row_id INT', 'url TEXT', 'missing_type TEXT', 'name TEXT'])
    
    for row in db.Execute('select * from nist_values'):
        skip_this_row = False
        row_comments = []

        row_id = row[0] # the first column is the ID, then fields are counted starting from 0
        row_dict = {}
        for title, colnum in titles2colnum.iteritems():
            row_dict[title] = row[1+colnum]
    
        # specific corrections to NIST
        if row_dict['URL'].find('&') == -1:
            continue
        url_id = row_dict['URL'].split('&')[1]
        if url_id == "T1=43KRE/EGG_1159": # NIST mistakenly list the pH as 1.4
            row_dict['pH'] = "7.4"
        if url_id == "T1=72WUR/HES_1276": # the alpha and beta appear in NIST but have been thrown away by our parsing
            row_dict['Reaction'] = "alpha-D-Glucose 6-phosphate(aq) = " + \
                "beta-D-Glucose 6-phosphate(aq)"
        if url_id == "T1=93VIN/GRU_1691": # NIST says nicotinamide mononucleotide instead of nicotinate mononucleotide
            row_dict['Reaction'] = "Nicotinate D-ribonucleotide(aq) + pyrophosphate(aq) = " + \
                "nicotinic acid(aq) + 5-Phospho-alpha-D-ribose 1-diphosphate(aq)"
        if url_id == "T1=73VEL/GUY_1167": # the concentration is actually in mM (not M)
            row_dict['c(MgCl2,mol dm-3)'] += " 10-3"
        if url_id == "T1=63GRE_1058":
            row_dict['Buffer'] = "potassium maleate (0.001 mol dm-3)" # originally NIST say 1.0 M, which is too high
        if url_id == "T1=69LAN/DEK_92":
            continue
        if url_id == "T1=98KIM/VOE_559":
            continue
        if url_id == "T1=74UEB/BLA_161" and row_dict['T(K)'] == '203.15': # typo in NIST, verified using original paper
            row_dict['T(K)'] = '303.15'
        if url_id == "T1=62GOL/WAG_1116": # Inconsistent with the original paper and with Alberty
            continue
        if url_id == "T1=62GOL/WAG_289": # Inconsistent with the original paper and with Alberty
            continue
        if url_id == "T1=80TER/RAB_994": # change ammonium carbamate to carbamate
            row_dict['Reaction'] = \
                "carbamate(aq) + H2O(l) = ammonia(aq) + carbon dioxide(aq)"
        if url_id == "T1=99KAT/UED_1558": # the value in the paper is different
            row_dict['K \''] = '4 10-4'
        if url_id == "T1=99KAT/UED_1557":
            row_dict['K'] = None
            if row_dict['pH'] == '9.0': # the value in the paper is different
                row_dict['K\''] = '4 10-3'
            if row_dict['pH'] == '7.4': # the value in the paper is different
                row_dict['K\''] = '8 10-3'
            if row_dict['pH'] == '6.0': # the value in the paper is different
                row_dict['K\''] = '3 10-2'
        if url_id == "T1=62HAL/FEN_1446": # replace (S)-methylmalonyl-CoA with methylmalonyl-CoA
            row_dict['Reaction'] = \
                "ATP(aq) + propanoyl-CoA(aq) + carbon dioxide(aq) = ADP(aq) + phosphate(aq) + methylmalonyl-CoA(aq)"
        if url_id == "T1=65KAZ/GRO_1447": # replace (S)-methylmalonyl-CoA with methylmalonyl-CoA
            row_dict['Reaction'] = \
                "ATP(aq) + propanoyl-CoA(aq) + carbon dioxide(aq) = ADP(aq) + phosphate(aq) + methylmalonyl-CoA(aq)"
        if url_id == "T1=59MIL/LUK_1225":
            row_dict['Reaction'] = "1-(5'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole" + \
                " = Fumarate + 1-(5'-Phosphoribosyl)-5-amino-4-imidazolecarboxamide"
        if url_id == "T1=88BED/HAD_1430": # reaction is not electron balanced (missing redox carried)
            continue
        if url_id == "T1=98LIA/QU_1582":
            row_dict['Reaction'] = "2 O2- + 2 H+ = O2 + H2O2"
        if url_id == "T1=61BER/BER_1431":
            row_dict['Reaction'] = "ATP + L-valine + tRNA(Val)" + \
                " = AMP + diphosphate + L-Valyl-tRNA(Val)"
        if url_id == "T1=48KOR_676":
            row_dict['Reaction'] = "ATP + beta-Nicotinamide mononucleotide = NAD + pyrophosphate"
        if url_id == "T1=89AIR_1032":
            row_dict['Reaction'] = "pantothenate + H2O = pantoic acid + beta-alanine"
        if url_id == "T1=58CAB/LEL_432":
            row_dict['Reaction'] = "UDP-glucose + D-glucose 6-phosphate = UDP + alpha,alpha'-trehalose 6-phosphate"
        if url_id == "T1=62HAY/NIS_525":
            row_dict['Reaction'] = "L-alanine + 3-oxopropanoate = beta-alanine + pyruvate" 
        if url_id == "T1=76GRE/BRI_328": # cytochrome c oxidase with O2, probably a NIST typo
            continue
        if url_id == "T1=66CAR/HUL_89" or url_id == "T1=66CAR/HUL_200": # the values for glycolate reductase conflict with 3 other papers
            continue
        if url_id == "T1=03RAN/VAN_1602": # In NIST, the reaction was written in the reverse direction, and mistakedly tagged this under K instead of K'
            row_dict["K'"] = u"0.001"
            row_dict["K"] = u""
        if url_id == "T1=82NG/WON_1590": # the value in the original paper is very different, and given as K'
            continue
        if url_id == "T1=95KAM/JUR_626": # value conflicts with two other papers and common sense (it is 10^9, which noone can measure)
            continue
        
        # NIST mistakedly tagged this under K instead of K' due to PDF rendering errors
        if url_id in ["T1=00BYR/GOL_1519", "T1=06XU/WES_1713",
                      "T1=07LIN/ALG_1584", "T1=00ROD/BAR_1603",
                      "T1=91PAR/HOR_1594", "T1=83HAA/KAR_1015",
                      "T1=83HAA/KAR_1013", "T1=79VAN_1688", "T1=99ELS_1527",
                      "T1=02TEW/HAW_1641"]:
            row_dict["K'"] = row_dict["K"]
            row_dict["K"] = u""
        
        # pH is missing
        if url_id == "T1=83HAA/KAR_1013":
            row_dict['pH'] = u"6.5"
            
        if url_id == "T1=06TAN/SUR_1609":
            # NIST made several mistakes here. First, they used K instead of K'.
            # Second, they miscalculated the value of K' by a factor of 10^12, since
            # it was given in uM^-1 and to convert it to M they divided by 10^6.
            # Third, they used the less reliable value of K' (0.05uM) instead of
            # 0.02 uM.
            row_dict["K'"] = u"5E-5"
            row_dict["K"] = u""
        
        new_row_dict = {}
        for old_title, new_title in title_mapping.iteritems():
            new_row_dict.setdefault(new_title, None)
            if row_dict[old_title] and not new_row_dict[new_title]:
                new_row_dict[new_title] = row_dict[old_title]
            elif row_dict[old_title] and new_row_dict[new_title]:
                raise Exception("Row %d (%s) in NIST has two values for %s" % 
                                (row_id, row_dict['URL'], new_title))

        new_row_dict['url'] = "http://xpdb.nist.gov/enzyme_thermodynamics/" + new_row_dict['url']
        new_row_dict['pMg'] = None
        
        if row_dict['cosolvent'] not in [None, 'none', '']:
            row_comments += [('cosolvent', row_dict['cosolvent'])]
            skip_this_row = True
        if row_dict['solvent'] not in [None, 'none', 'H2O', '']:
            row_comments += [('solvent', row_dict['solvent'])]
            skip_this_row = True
        
        # specific corrections to NIST
        new_row_dict['reaction'] = new_row_dict['reaction'].replace(' +-D-', ' + D-')
        new_row_dict['reaction'] = new_row_dict['reaction'].replace(' +-lipoate', ' + lipoate')
        try:
            new_row_dict['kegg_reaction'] = parse_reaction(new_row_dict['reaction'],
                                                       compound_aliases)
        except MissingCompoundsFromKeggException as e:
            for name in e.names:
                row_comments += [('reactant', name)]
            skip_this_row = True

        for key in real_titles:
            try:
                new_row_dict[key] = make_sure_it_is_float(new_row_dict[key])
            except ValueError:
                raise Exception("Cannot parse row %d, the value of %s is %s" %
                                 (row_id, key, new_row_dict[key]))
        
        ### Salts ###
        Mg_conc = 0
        charge_conc_pairs = []
        for title, cc_lists in salt_titles.iteritems():
            try:
                conc = make_sure_it_is_float(row_dict[title])
            except ValueError:
                continue
            if not conc:
                continue
            
            if title in ['pMn', 'pMg']:
                conc = 10**(-conc)
                
            if title in Mg_titles:
                Mg_conc += conc

            if type(cc_lists) == types.StringType:
                base_charge, pKa_list = buffer_dict[cc_lists]
                charge_conc_pairs += get_buffer_charges(base_charge, pKa_list, 
                                                        conc, new_row_dict['pH'])
            else:
                if not cc_lists:
                    raise Exception(title)
                for charge, coeff in cc_lists:
                    charge_conc_pairs.append((charge, conc*coeff))
                    
        ### Buffers ###
        reaction_buffer = choose_buffer(row_dict)
        missing_buffers = {}
        charge_conc_pairs += buffer_match(reaction_buffer, 
            buffer_dict, new_row_dict['pH'], missing_buffers)
        if missing_buffers:
            row_comments += [('buffers', str(missing_buffers))]
        
        if new_row_dict['I'] is None and charge_conc_pairs:
            # on 06/05/2012 Elad & Arren decided to stop using this calculation since 
            # it ignores the charges of the reactants themselves, and thus 
            # usually underestimates I by a large margin.
            #
            #new_row_dict['I'] = np.round(0.5 * sum([(conc * ch**2) for (ch, conc) in charge_conc_pairs]), 2)
            pass
        
        if Mg_conc > 0:
            new_row_dict['pMg'] = -np.round(np.log10(Mg_conc), 2)

        if skip_this_row:
            for missing_type, name in row_comments:
                db.Insert('nist_errors', [row_id, new_row_dict['url'], missing_type, name])
        else:
            db.Insert('nist_equilibrium', [new_row_dict[k] for k in new_titles])

    for row in csv.DictReader(open(NIST_ADDITIONS_FILE, 'r')):
        for key in row.keys():
            if row[key] == '':
                row[key] = None
        db.Insert('nist_equilibrium', [row[k] for k in new_titles])

    db.Commit()

    for row in db.Execute("SELECT COUNT(*) FROM nist_equilibrium "
                          "WHERE K_tag IS NOT NULL OR K IS NOT NULL"):
        print "%d reactions have been succesfully mapped!" % row[0]

if __name__ == "__main__":
    db = SqliteDatabase('../data/public_data.sqlite')
    WriteDataToDB(db)
    
    #test_buffer_methods()