import re
import csv
from math import log10
import pylab
from pygibbs.thermodynamic_constants import default_pH

BUFFERS_CSV_FNAME = '../data/thermodynamics/pKa_of_buffers.csv'
NIST_CSV_FNAME = '../data/thermodynamics/nist_equilibrium_raw.csv'
KEGG_COMPOUND_FNAME = '../data/thermodynamics/nist_compounds_to_kegg.csv'
OUTPUT_FNAME = '../data/thermodynamics/nist.csv'

class MissingCompoundsFromKeggException(Exception):
    def __init__(self, names):
        self.names = names
    def __str__(self):
        return "Could not find this reactant in the list: " + \
            ', '.join(self.names)

def remove_superfluous_chars(s):
    to_remove = ['(g)', '(l)', '(aq)', '(liq)', '(s)', '(sln)', '"']
    for substr in to_remove:
        s = s.replace(substr, '')
    s = s.strip()
    if s[0] == '-':
        s = s[1:]
    return s.lower()

def parse_single_reactant(s, compound_aliases):
    s = remove_superfluous_chars(s)
    if s in compound_aliases:
        if s == "carbon dioxide":
            # change carbon dioxide changed to HCO3
            return "C00288"
        else:
            return compound_aliases[s]

    tmp = re.findall('^(\d+) (.*)', s)
    if tmp:
        count, name = tmp[0]
        if name in compound_aliases:
            return count + " " + compound_aliases[name]
    
    raise MissingCompoundsFromKeggException([s])                

def parse_reaction_side(s, compound_aliases):
    res = []
    missing_names = []
    for r in s.split(' + '):
        try:
            res.append(parse_single_reactant(r, compound_aliases))
        except MissingCompoundsFromKeggException as e:
            missing_names += e.names
    if missing_names:
        raise MissingCompoundsFromKeggException(missing_names)
    else:
        return " + ".join(res)

def parse_reaction(s, compound_aliases):
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
    if not s or s == 'not given':
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
        if re.findall('\d+\.\d+', row_dict[key]):
            return row_dict[key].lower()
    return None

def get_buffer_charges(base_charge, pKa_list, conc, pH):
    if not pH:
        pH = default_pH
    
    pKa_list = sorted(pKa_list, reverse=True)
    
    species_proportions = []
    for n in xrange(len(pKa_list) + 1):
        p = pylab.prod([10**(pKa - pH) for pKa in pKa_list[:n]])
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

def get_buffer_charges2(base_charge, pKa_list, conc, pH):
    i = 0
    charge = base_charge
    while i < len(pKa_list) and pH > pKa_list[i]:
        i += 1
        charge -= 1
    return [(charge, conc)]

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
    for row in csv.reader(open(KEGG_COMPOUND_FNAME, 'r')):
        cid = int(row[0])
        for alias in row[2].split(';'):
            alias = remove_superfluous_chars(alias)
            compound_aliases[alias] = "C%05d" % cid
    
    #manually added compounds
    compound_aliases["1-dodecanoic acid"] = "C02679"
    compound_aliases[",-trehalose 6-phosphate"] = "C00689"
    compound_aliases["nad(ox)"] = "C00003"
    compound_aliases["nad(red)"] = "C00004"
    compound_aliases["nadp(ox)"] = "C00006"
    compound_aliases["1/2 o2"] = "0.5 C00007"
    compound_aliases["gdpglucose"] = "C00394"
    compound_aliases["adpglucose"] = "C00498"
    compound_aliases["dadp"] = "C00206"
    compound_aliases["3-oxobutanoate"] = "C00109"
    compound_aliases["pyrophosphate  this is an approximate result."] = "C00013"
    compound_aliases["vitamin a alcohol"] = "C00473"
    compound_aliases["n-decanoic acid"] = "C01571"
    compound_aliases["2'-deoxyinosine"] = "C05512"
    compound_aliases["n-octonic acid"] = "C06423"
    compound_aliases["orotidine 5 '-phosphate"] = "C01103"
    compound_aliases["5-pregnane-17,21-diol-3,11,20-trione"] = "C05469"
    compound_aliases["5-pregnane-3,17,21-triol-11,20-dione"] = "C05471"
    compound_aliases["5-phospho--d-ribose 1-diphosphate"] = "C00119"
    compound_aliases["nadp(red)"] = "C00005"
    compound_aliases["n-w-phospho-l-arginine"] = "C05945"
    compound_aliases["6-amino-d-glucose 6-phosphate"] = "C00352"
    compound_aliases["cis-cis-hexadienedioate"] = "C02480"
    compound_aliases["n-octanoic acid"] = "C06423"
    compound_aliases["glycinamide"] = "C03838"
    compound_aliases["lyxose"] = "C01508"
    compound_aliases["adenosine 3':5'-(cyclic)phosphate"] = "C00575"
    compound_aliases["l-arginyl-trna-arg"] = "C02163"
    compound_aliases["/-)-1-phenyl-1-ethanol"] = "C07112"
    compound_aliases["udpgalactose"] = "C00052"
    compound_aliases["d-1-piperidine-2-carboxylate"] = "C04092"
    compound_aliases["gdpmannose"] = "C00096"
    compound_aliases["benzyloxycarbonylglycine"] = "C03710"  
    compound_aliases["l-tyrosyl-trna-tyr"] = "C02839"
    compound_aliases["guanosine 5 '-phosphate"] = "C00144"
    compound_aliases["nomega-phospho-l-arginine"] = "C05945"
    compound_aliases["adp--s"] = "C01469"
    compound_aliases["\)-1-phenyl-1-propyl acetate"] = "C17666"
    compound_aliases["adenylyl-[l-glutamate:ammonia ligase(adp-forming)]"] = "C01299"
    compound_aliases[")-1-phenyl ethanol"] = "C07112"
    compound_aliases["dihydro--lipoate"] = "C02147"
    compound_aliases["trna-tyr"] = "C00787"
    compound_aliases["phosphotaurocyamine"] = "C03149"
    compound_aliases["trna-arg"] = "C01636"
    compound_aliases["d-glucose 6-phosphate-2-"] = "C00092"
    compound_aliases["1-phenyl-1-ethanone"] = "C07113"
    compound_aliases[",-trehalose"] = "C01083"
    compound_aliases["(r)-3-phosphoglycerate"] = "C00197"
    compound_aliases["[l-glutamate:ammonia ligase(adp-forming)]"] = "C01281"
    compound_aliases["1,2,3-trioctanoyl glycerol"] = "C013044" #trioctanoylglycerol
    compound_aliases["trna-phe"] = "C01648"
    compound_aliases["trna-lys"] = "C01646"
    compound_aliases["trna-his"] = "C01643"
    compound_aliases["trna-ser"] = "C01650"
    compound_aliases["trna-ile"] = "C01644"
    compound_aliases["l-isoleusine"] = "C00407"
    compound_aliases["(n-1) h2o"] = "(n-1) C00001"
    compound_aliases["estradiol-17"] = "C00951"
    compound_aliases["trna-thr"] = "C01651"
    compound_aliases["trna-val"] = "C01653"
    compound_aliases["h2o(1)"] = "C00001"
    compound_aliases["l-arginosuccinate"] = "C03406"
    compound_aliases["glucose 1-phosphate"] = "C00103"
    compound_aliases["pyruvate-"] = "C00022"
    compound_aliases["l-lysyl-trna-lys"] = "C01931"
    compound_aliases["l-histidyl-trna-his"] = "C02988"
    compound_aliases["l-seryl-trna-ser"] = "C02553"
    compound_aliases["l-isoleucyl-trna-ile"] = "C03127"
    compound_aliases["5-amino-1--d-ribosyl-4-imidazolecarboxamide 5'-phosphate"] = "C04677"
    compound_aliases["3'-phosphoadenylylsulfate"] = "C00053"
    compound_aliases["5-androstane-3,17-dione"] = "C00674"
    compound_aliases["L--glycerophosphate(aq)"] = "C00093"
    compound_aliases["glutathione (reduced)"] = "C00051"
    compound_aliases["glutathoine(red)"] = "C00051"
    compound_aliases["dimethylmaleate"] = "C00922"
    compound_aliases["(s)-proline"] = "C00148"
    compound_aliases["octanal"] = "C01545"
    compound_aliases["orthophospate"] = "C00009"
    compound_aliases["d-ribulose 1,5-biphosphate"] = "C01182"
    compound_aliases["heptanoate"] = "C17714"
    compound_aliases["propanonyl-coa"] = "C00100"
    compound_aliases["(r)-glyceraldehyde"] = "C00577"
    compound_aliases["ammonium carbamate"] = "C01563" # AKA carbamate
    compound_aliases["tributyrylglycerol"] = "C13870" # AKA Glyceryl tributyrate
    compound_aliases["3-oxobutanoyl-coa"] = "C00332" # AKA acetoacetyl-CoA
    compound_aliases["ap-nadh"] = "C90001" # 3-acetyl NADH
    compound_aliases["ap-nad"] = "C90002" # 3-acetyl NAD+
    
    return compound_aliases

# rows with K instead of K'
reverse_transformed_rows = \
['T1=00BYR/GOL_1519',
'T1=00ROD/BAR_1603',
'T1=01TIS/IHL_1682',
'T1=02TEW/HAW_1641',
'T1=03RAN/VAN_1602',
'T1=06TAN/SUR_1609',
'T1=06XU/WES_1713',
'T1=07LIN/ALG_1584',
'T1=76CRA/WAI_1521',
'T1=79VAN_1688',
'T1=82NG/WON_1590',
'T1=91PAR/HOR_1594',
'T1=98DIE/STR_1524',
'T1=99ELS_1527',
'T1=99KAT/UED_1557',
'T1=99NIE/SCH_1592']

################################################################################

def rewrite_nist(buffer_dict, compound_aliases):
    """
        each reaction is composed of 'substrates' and 'products'.
        each compound from both classes is matched to aliases and recieve its 
        representative name
    """
    successful_row_counter = 0
    missing_compounds = {}
    missing_buffers = {}

    nist_reader = csv.reader(open(NIST_CSV_FNAME, 'r'))
    titles = nist_reader.next()
    salt_titles = []
    for t in titles[14:46]:
        tmp = re.findall('[cm]\(([\w\d\+\-]+),mol[\s\.][dk][mg]-[13]\)', t)
        if tmp:
            salt_titles.append(tmp[0])
            continue
        tmp = re.findall('[cm]\(([\w\d\+\-]+)\)', t)
        if tmp:
            salt_titles.append(tmp[0])
            continue
        salt_titles.append(t)
    
    new_titles = ['URL', 'Reference_id', 'Method', 
        'Evaluation', 'EC value', 'Enzyme', 'KEGG Reaction', 'Reaction',
        'K', 'Temp', 'Ionic strength', 'pH', 'pMg', 'comment', 'Ktype']
    
    nist_writer = csv.writer(open(OUTPUT_FNAME, 'w'))
    nist_writer.writerow(['url', 'reference_id', 'method', 'evaluation', 'ec',
        'enzyme', 'kegg_reaction', 'reaction', 'K', 'T', 'I', 'pH', 'pMg', 
        'comment', 'Ktype'])
    
    row_counter = 0
    for row in nist_reader:
        row_counter += 1
    
        # specific corrections to NIST
        if row[0].find('&') == -1:
            continue
        url_id = row[0].split('&')[1]
        if url_id == "T1=43KRE/EGG_1159": # NIST mistakenly list the pH as 1.4
            row[titles.index('pH')] = "7.4"
        if url_id == "T1=72WUR/HES_1276": # the alpha and beta appear in NIST but have been thrown away by our parsing
            row[titles.index('Reaction')] = "alpha-D-Glucose 6-phosphate(aq) = " + \
                "beta-D-Glucose 6-phosphate(aq)"
        if url_id == "T1=93VIN/GRU_1691": # NIST says nicotinamide mononucleotide instead of nicotinate mononucleotide
            row[titles.index('Reaction')] = "Nicotinate D-ribonucleotide(aq) + pyrophosphate(aq) = " + \
                "nicotinic acid(aq) + 5-Phospho-alpha-D-ribose 1-diphosphate(aq)"
        if url_id == "T1=73VEL/GUY_1167": # the concentration is actually in mM (not M)
            row[titles.index('c(MgCl2,mol dm-3)')] += " 10-3"
        if url_id == "T1=63GRE_1058":
            row[titles.index('Buffer')] = "potassium maleate (0.001 mol dm-3)" # originally NIST say 1.0 M, which is too high
        if url_id == "T1=69LAN/DEK_92":
            continue
        if url_id == "T1=98KIM/VOE_559":
            continue
        if url_id == "T1=74UEB/BLA_161" and row[titles.index('Temp')] == '203.15': # typo in NIST, verified using original paper
            row[titles.index('Temp')] = '303.15'
        if url_id == "T1=62GOL/WAG_1116": # possible huge mistake in NIST
            continue
        if url_id == "T1=80TER/RAB_994": # type on NIST, where it's written 2 ammonia instead of only 1
            row[titles.index('Reaction')] = \
                "ammonium carbamate(aq) + H2O(l) = ammonia(aq) + carbon dioxide(aq)"
        
        row_dict = dict([(titles[i], row[i]) for i in xrange(len(titles))])
        row_dict['URL'] = "http://xpdb.nist.gov/enzyme_thermodynamics/" + row_dict['URL']
        row_dict['KEGG Reaction'] = None
        row_dict['comment'] = None
        row_dict['K'] = make_sure_it_is_float(row_dict['K'])
        if url_id in reverse_transformed_rows:
            row_dict["Ktype"] = "K"
        else:
            row_dict["Ktype"] = "K\'"
    
        if row_dict['Reaction'].find('=') != -1:
            # specific corrections to NIST
            row_dict['Reaction'] = row_dict['Reaction'].replace(' +-D-', ' + D-')
            row_dict['Reaction'] = row_dict['Reaction'].replace(' +-lipoate', ' + lipoate')
            try:
                row_dict['KEGG Reaction'] = parse_reaction(row_dict['Reaction'],
                                                           compound_aliases)
            except MissingCompoundsFromKeggException as e:
                if row_dict['K']: # count missing compounds only if the reaction has an observed Keq
                    for s in e.names:
                        missing_compounds[s] = missing_compounds.get(s, 0) + 1
        
        for key in ['Temp', 'Ionic strength', 'pH']:
            row_dict[key] = make_sure_it_is_float(row_dict[key])
    
        if row_dict['x(cosolvent)']:
            row_dict['comment'] = 'cosolvent => ' + row_dict['x(cosolvent)']
        if row_dict['cosolvent']:
            row_dict['comment'] = 'cosolvent => ' + row_dict['cosolvent']
    
        Mg_conc = 0
        charge_conc_pairs = []
        ### Salts ###
        concentrations = row[14:46]
        for i in xrange(len(concentrations)):
            try:
                conc = make_sure_it_is_float(concentrations[i])
            except ValueError:
                continue
            if not conc:
                continue
            if salt_titles[i] in ['KCl', 'NaCl']:
                charge_conc_pairs += [(1, conc), (-1, conc)]
            elif salt_titles[i] == 'Na+':
                charge_conc_pairs += [(1, conc)]
            elif salt_titles[i] in ['Mg2+', 'Mg']:
                charge_conc_pairs += [(2, conc)]
            elif salt_titles[i] in ['MgCl2', 'MnCl2', 'ZnCl2', 'CaCl2']:
                charge_conc_pairs += [(2, conc), (-1, 2*conc)]
            elif salt_titles[i] in ['MnSO4', 'MgSO4']:
                charge_conc_pairs += [(2, conc), (-2, conc)]
            elif salt_titles[i] in ['pMn', 'pMg']:
                charge_conc_pairs += [(2, 10**(-conc))]
            elif salt_titles[i] in ['phosphate', 'orthophosphate']:
                base_charge, pKa_list = buffer_dict['phosphate']
                charge_conc_pairs += get_buffer_charges(base_charge, pKa_list, 
                                                        conc, row_dict['pH'])
            elif salt_titles[i] in ['Tris']:
                base_charge, pKa_list = buffer_dict['tris']
                charge_conc_pairs += get_buffer_charges(base_charge, pKa_list, 
                                                        conc, row_dict['pH'])
                
            if salt_titles[i] in ['Mg2+', 'Mg', 'MgCl2', 'MgSO4']:
                Mg_conc += conc
            elif salt_titles[i] == 'pMg':
                Mg_conc += 10**(-conc)
    
        ### Buffers ###
        row_dict['buffer'] = choose_buffer(row_dict)
        charge_conc_pairs += buffer_match(row_dict['buffer'], 
                                          buffer_dict, row_dict['pH'], missing_buffers)
        
        if row_dict['pH']:
            row_dict['pH'] = '%.3g' % row_dict['pH']
        if not row_dict['Ionic strength']:
            row_dict['Ionic strength'] = "%.3g" % (0.5 * sum([(conc * ch**2) for (ch, conc) in charge_conc_pairs]))
        if Mg_conc > 0:
            row_dict['pMg'] = "%.3g" % -log10(Mg_conc)
            
        nist_writer.writerow([row_dict[k] for k in new_titles])
        if row_dict['K'] and row_dict['KEGG Reaction'] and not row_dict['comment']:
            successful_row_counter += 1
    return missing_compounds, missing_buffers, successful_row_counter

def main():
    buffer_dict = load_buffer_dict()
    compound_aliases = load_compound_aliases()
    missing_compounds, missing_buffers, successful_row_counter = rewrite_nist(buffer_dict, compound_aliases)
    
    for name, count in sorted(missing_buffers.iteritems(), key=lambda x:x[1]):
        print "Missing buffer: %s (%d times)" % (name, count) 

    for name, count in sorted(missing_compounds.iteritems(), key=lambda x:x[1]):
        print "Missing compound: %s (%d times)" % (name, count)
        
    print "Deciphered %d reactions!" % successful_row_counter

def test_buffer_methods():
    pH_range = pylab.arange(1, 14.01, 0.1)
    N = len(pH_range)
    base_charge = 10
    pKa_list = [5, 8, 11]
    c_total = 0.1
    charge = pylab.zeros((N, 2))
    for i in xrange(N):
        charge[i, 0] = 0.5 * sum([(ch**2 * conc) for (ch, conc) in get_buffer_charges(base_charge, pKa_list, c_total, pH_range[i])])
        charge[i, 1] = 0.5 * sum([(ch**2 * conc) for (ch, conc) in get_buffer_charges2(base_charge, pKa_list, c_total, pH_range[i])])
    pylab.plot(pH_range, charge)
    pylab.xlabel('pH')
    pylab.ylabel('Ionic strength')
    pylab.show()

if __name__ == "__main__":
    main()
    #test_buffer_methods()