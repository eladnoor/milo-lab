from scipy import stats
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs import thermodynamics
from pygibbs.thermodynamic_constants import R, default_I, default_pH,\
    default_pMg, default_T
from pygibbs.groups import GroupContribution
import pylab
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggNonCompoundException,\
    KeggReactionNotBalancedException
from pygibbs.kegg_reaction import Reaction
from toolbox import plotting
from SOAPpy import WSDL 
from pygibbs.metacyc import MetaCyc, MetaCycNonCompoundException
from toolbox import util
import logging
import matplotlib
import re
import numpy
import random
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy

def try_kegg_api():
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/dG0_test.html')
    G = GroupContribution(db, html_writer=html_writer)
    G.init()
    
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    
    rid_file = open('../res/eco_rids.txt', 'w')
    rids = set()
    for x in serv.list_pathways('eco'):
        pathway_id = x['entry_id']
        for reaction_id in serv.get_reactions_by_pathway(pathway_id):
            rid = int(reaction_id[4:])
            if rid not in rids:
                rids.add(rid)
                rid_file.write('%d\n' % rid)
    rid_file.close()
            
    c_mid = 1e-3
    pH, pMg, I, T = (7.0, 3.0, 0.1, 298.15)
    
    rid2reversibility = {}
    misses = 0
    for rid in sorted(rids):
        try:
            sparse = G.kegg.rid2sparse_reaction(rid)
            r = CalculateReversability(sparse, G, c_mid, pH, pMg, I, T)
            rid2reversibility[rid] = r
        except thermodynamics.MissingCompoundFormationEnergy:
            misses += 1
            continue
    
    print 'hits = %d, misses = %d' % len(rid2reversibility), misses
    median = pylab.median(rid2reversibility.values())
    print 'median = %.1f' % median

    pylab.figure()
    pylab.hold(True)
    plotting.cdf(rid2reversibility.values(), 'all reactions', 'r', show_median=True)
    pylab.show()
    

WATER = 1
HPLUS = 80


def GetConcentrationMap():
    kegg = Kegg.getInstance()
    cmap = {}    
    for cid in kegg.get_all_cids():
        lower, upper = kegg.get_bounds(cid)
        if lower and upper:
            # In the file we got this data from lower = upper 
            cmap[cid] = lower
    cmap[WATER] = 1
    cmap[HPLUS] = 1
    return cmap



def GetFullConcentrationMap(G):
    G.read_compound_abundance('../data/thermodynamics/compound_abundance.csv')
    cmap = {}
    kegg_handle = Kegg.getInstance()
    
    for cid in kegg_handle.get_all_cids():
        cmap[cid] = G.get_concentration(cid, media='glucose', c0=None)
        if cmap[cid] == None:
            lower, upper = kegg_handle.get_bounds(cid)
            if lower and upper:
                # In the file we got this data from lower = upper 
                cmap[cid] = lower
            else:
                cmap.pop(cid)
    cmap[WATER] = 1
    cmap[HPLUS] = 1
    return cmap

def ConcentrationFactor(sparse_reaction,
                        concentration_map,
                        c_mid):
    factor = 0.0
    for cid, stoic in sparse_reaction.iteritems():
        concentration = concentration_map.get(cid, None) or c_mid
        factor += pylab.log(concentration) * stoic
    return factor

def CalculateReversability(sparse, thermo, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T,
                           concentration_map=None):
    cmap = concentration_map or {}
    rxn = Reaction("Unknown", sparse)
    dG0 = thermo.reaction_to_dG0(rxn, pH, pMg, I, T)
    
    ln_Gamma = ConcentrationFactor(sparse, cmap, c_mid)
    sum_abs_s = sum([abs(x) for k, x in sparse.iteritems()
                     if k not in cmap])
    if sum_abs_s == 0:
        return None
    else:
        return -(2 / (sum_abs_s * pylab.log(10))) * (ln_Gamma + dG0/(R*T))
        

def CalculateReversabilityV2(sparse, thermo, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T,
                           concentration_map=None):
    cmap = concentration_map or {}
    dG0 = thermo.reaction_to_dG0(sparse, pH, pMg, I, T)
    
    cfactor = ConcentrationFactor(sparse, cmap, c_mid)
    
    if (cfactor == 0):
        return None
    else:
        return  (-dG0/(R * T * cfactor)) / pylab.log(10)

def CalculateReversabilitydeltaG(sparse, thermo, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T,
                           concentration_map=None):
    cmap = concentration_map or {}
    dG0 = thermo.reaction_to_dG0(sparse, pH, pMg, I, T)
    
    cfactor = ConcentrationFactor(sparse, cmap, c_mid)
    
    return dG0 + R * T * cfactor

def calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, cmap, id):
    kegg = Kegg.getInstance()
    histogram = {}
    histogram['Not first'] = []
    histogram['Rest'] = []
    
    rel_histogram = {}
    rel_histogram['Not first'] = []
    
    n_first_max = 0
    n_valid = 0
    
    hits = 0
    misses = 0
    total_rxns = 0
    only_currency = 0
    non_balanced = 0
    
    n_pathways_with_thermo_data = 0
    n_short_pathways = 0
    
    rxns_map = {}
    total_rxns_map = {}
    rxn_reversibilities = []
    
    debug_file = open('../res/kegg_' + id + '_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Module\tPosition\tReaction Name\tDefinition\tEC list\tEquation\tRev IND\n")
    
    for mid, rid_flux_list in kegg.mid2rid_map.iteritems():
        if not rid_flux_list or len(rid_flux_list) < 2:
            n_short_pathways += 1
            continue
        
        pw_r_map = {}
        pw_r_sum = 0
        pw_n_r = 0
        
        has_thermo_data = 0
        for i, (rid, flux) in enumerate(rid_flux_list):
            if rid not in total_rxns_map:
                total_rxns_map[rid] = 1
                total_rxns += 1
                
            try:
                reaction = kegg.rid2reaction(rid)
                reaction.Balance(balance_water=True)
                sparse = reaction.sparse
                
                # delta G r = CalculateReversabilitydeltaG(sparse, G, c_mid, pH, pMg, I, T, concentration_map=cmap)
                r = CalculateReversability(sparse, G, c_mid, pH, pMg, I, T,
                                           concentration_map=cmap)
                if r == None:
                    if rid not in rxns_map:
                        rxns_map[rid] = 1
                        only_currency += 1
                    continue
                
                r *= flux
                
                rxn_reversibilities.append((abs(r), rid))
                has_thermo_data = 1
                #r = abs(r)
                
                rxn = kegg.rid2reaction_map[rid]
                dbg = "%s__DEL__%d__DEL__%s__DEL__%s__DEL__%s__DEL__%s__DEL__%f\n" % (kegg.mid2name_map[mid], i+1, rxn.name, rxn.definition, rxn.ec_list, rxn.equation, r)
                dbg = re.sub ('\t', ' ', dbg)
                debug_file.write(re.sub('__DEL__', '\t', dbg))
                
                if i > 0:
                    histogram['Not first'].append(r)
                if i > 4:
                    histogram['Rest'].append(r)
                else:
                    histogram.setdefault(i+1, []).append(r)
                
                if rid not in rxns_map:
                    rxns_map[rid] = 1
                    hits += 1
                
                pw_n_r += 1
                pw_r_sum += abs(r)
                if (i in pw_r_map):
                    pw_r_map[i] = pw_r_map[i] + ';' + str(r)
                else:
                    pw_r_map[i] = str(r)
                            
            except thermodynamics.MissingCompoundFormationEnergy:
                if rid not in rxns_map:
                    rxns_map[rid] = 1
                    misses += 1
                continue
            except KeggReactionNotBalancedException:
                    if rid not in rxns_map:
                        rxns_map[rid] = 1
                        non_balanced += 1
                    #print 'Reaction cannot be balanced, uid: %s' % rxn
                    continue
        
        if has_thermo_data == 1:
            n_pathways_with_thermo_data += 1
                
        if (pw_n_r >= 1):
            avg_abs_r = pw_r_sum / pw_n_r
            if (len(pw_r_map) > 1 and 0 in pw_r_map):
                first_r = -1000
                curr_max = -1000
                n_valid += 1
                for pos,r_list in pw_r_map.iteritems():
                    for r in r_list.split(';'):
                        r = float(r)
                        if (pos == 0):
                            first_r = r
                        curr_max = r if r > curr_max else curr_max
                        r = float(r) / avg_abs_r
                        rel_histogram.setdefault(pos+1, []).append(r)
                        if pos > 0:
                            rel_histogram['Not first'].append(r)
                if (first_r == curr_max):
                        n_first_max += 1
    
    rxn_reversibilities.sort(reverse=True)
    print rxn_reversibilities[:20]
                                        
    debug_file.close()
    
    debug_file = open('../res/kegg_' + id + '_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev_stats.txt', 'w')
    debug_file.write("Total number of pathways with positional data: %d\n" % len(kegg.mid2rid_map))
    debug_file.write("Single reaction or pathways without reactions data: %d\n" % n_short_pathways)
    debug_file.write("Total number of pathways with thermodynamic data: %d\n" % n_pathways_with_thermo_data)
    
    debug_file.write("Total unique reactions: %d\n" % total_rxns)
    debug_file.write("Reactions with known dG0: %d\n" % hits)
    debug_file.write("Reactions with unknown dG0: %d\n" % misses)
    debug_file.write("Non balanced reactions: %d\n" % non_balanced)
    debug_file.write("Reactions without unknown concentrations: %d\n" % only_currency)
    debug_file.close()
    
    return histogram,rel_histogram,(float(n_first_max)/n_valid)

def plot_histogram(histogram, html_writer, title='', max_pathway_length=8, xmin=None, xlim=20, error_bars=True, min_to_show=20, legend_loc='upper left'):
    fig = pylab.figure()

    pylab.hold(True)

    reps = 1000
    
    y_offset = 0
    offset_step = 0.007
    colors = {1:'r', 2:'orange', 3:'green', 4:'cyan', 5:'blue', 'Rest':'violet', 'Not first':'k--', 'No known regulation':'grey', 'Activated':'green', 'Inhibited':'r', 'Mixed regulation':'blue'}
    for key, value in histogram.iteritems():
        if len(value) >= min_to_show:
            m = stats.cmedian(value)
            
            sample_std = None
            
            if error_bars:
                sample_vals = []
                i = 0
                while i < reps:
                    samples = []
                    while len(samples) < len(value):
                        samples.append(random.choice(value))
                    sample_vals.append(pylab.median(samples))
                    i += 1
                
                sample_std = pylab.std(sample_vals)
                        
            plotting.cdf(value, label='%s (med=%.1f, N=%d)' % \
                (key, m, len(value)),
                style=colors.get(key, 'grey'), std=sample_std, y_offset=y_offset)
            y_offset += offset_step
            

    xmin = -1 * xlim if xmin == None else xmin
    pylab.xlim(xmin, xlim)
    pylab.xlabel('Irreversability')
    #pylab.xlabel('deltaG')
    pylab.ylabel('Cumulative distribution')
    legendfont = matplotlib.font_manager.FontProperties(size=11)
    pylab.legend(loc=legend_loc, prop=legendfont)
    pylab.title(title)
    pylab.hold(False)
    
    if 'Not first' in histogram:
        print '%s, first vs. non-first ranksum test: ' % title + '(%f, %f)' % stats.ranksums(histogram[1], histogram['Not first'])
    
    if 'Inhibited' in histogram:
        print '%s, inhibited vs. non-regulated ranksum test: ' % title + '(%f, %f)' % stats.ranksums(histogram['Inhibited'], histogram['No known regulation'])
         
    
    #for k1, h1 in histogram.iteritems():
    #    for k2, h2 in histogram.iteritems():
    #        print k1, k2, stats.ranksums(h1, h2)
    
    return fig

def plot_bars (pos_count, title='', max_pathway_length=8, legend_loc='upper right'):
    n_labels = len(pos_count)
    ind = numpy.arange(max_pathway_length)
    width = 0.2
    
    fig = pylab.figure()
    pylab.hold(True)
    ax = fig.add_subplot(111)

    colors = {'No known regulation':'grey', 'Activated':'green', 'Inhibited':'red', 'Mixed regulation':'blue'}
    plot_order = ['Inhibited', 'Mixed regulation', 'Activated', 'No known regulation']    
    
    i = 0
    for label in plot_order:
        curr_vals = pos_count[label][1:max_pathway_length+1]
        if (sum(curr_vals) < 20):
            n_labels -= 1
            continue
        ax.bar(ind + i * width, tuple([j * 1.0 /sum(curr_vals) for j in curr_vals]), width, color=colors[label], label=('%s (%d)' % (label, sum(curr_vals))))
        i += 1
    
    ax.set_ylabel('Fraction of reactions per type')
    ax.set_xlabel('Position in pathway')
    
    ax.set_xticks(ind+ width * n_labels/2)
    ax.set_xticklabels( ind + 1 )
    
    legendfont = matplotlib.font_manager.FontProperties(size=11)
    pylab.legend(loc=legend_loc, prop=legendfont)
    pylab.title(title)

    pylab.hold(False)
    
    return fig

def plot_stacked_bars (pos_count, title='', max_pathway_length=8, legend_loc='upper right'):
    
    ind = numpy.arange(max_pathway_length)
    width = 0.35
    
    fig = pylab.figure()
    pylab.hold(True)
    ax = fig.add_subplot(111)

    colors = {'No known regulation':'grey', 'Activated':'green', 'Inhibited':'red', 'Mixed regulation':'blue'}
    plot_order = ['Inhibited', 'Mixed regulation', 'Activated', 'No known regulation']
    
    total = numpy.array([0] * max_pathway_length)
    for label,vals in pos_count.iteritems():
        total += numpy.array(vals[1:max_pathway_length+1])

    so_far = numpy.array([0] * max_pathway_length)
    for label in plot_order:
        vals = pos_count[label]
        curr_vals = numpy.array(vals[1:max_pathway_length+1]) * 1.0 / total
        ax.bar(ind, tuple(curr_vals), width, color=colors[label], bottom=tuple(so_far), label=('%s (%d)' % (label, sum(vals[1:max_pathway_length+1]))))
        so_far = so_far + numpy.array(curr_vals)
    
    ax.set_ylabel('Perc. of type per position')
    ax.set_xlabel('Position in pathway')
    
    ax.set_xticks(ind + 0.5 * width)
    ax.set_xticklabels( ['%d (%d)' % (ind[i] + 1, total[i]) for i in ind] )
    
    legendfont = matplotlib.font_manager.FontProperties(size=11)
    pylab.legend(loc=legend_loc, prop=legendfont)
    pylab.title(title)

    pylab.hold(False)
    
    return fig
    
def plot_bootstrap_stats(histogram, title=''):
    fig = pylab.figure()
    
    pylab.hold(True)
    
    colors = {1:'r', 2:'orange', 3:'green', 4:'cyan', 5:'blue', 'Rest':'violet', 'Not first':'k'}
    
    plotting.bootstrap(histogram, colors)
        
    pylab.title(title)
    
    pylab.hold(False)
    
    return fig



def calculate_metacyc_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, metacyc, cmap, id):
    
    kegg = Kegg.getInstance()
    histogram = {}
    histogram['Not first'] = []
    histogram['Rest'] = []

    rel_histogram = {}
    rel_histogram['Not first'] = []
    
    reg_hist = {}
    
    n_first_max = 0
    n_valid = 0
    
    hits = 0
    misses = 0
    comp_translation_misses = 0
    flux_misses = 0
    rxn_map_misses = 0
    meta_compound_misses = 0
    meta_non_balanced = 0
    n_only_currency = 0
    
    n_pathways_with_thermo_data = 0
    n_short_pw = 0    
    
    rxns_map = {}
    total_rxns_map = {}
    
    total_rxns = 0
    
    debug_file = open('../res/metacyc_' + id + ('_constrained_' if len(cmap) > 0 else '_non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Pathway\tPosition\tReaction Name\tEC list\tEquation\tRev IND\n")
    
    for pathway in metacyc.uid2pathway_map.itervalues():
        rxns_dict = pathway.GetRxnsOrder()
        rxn_dirs = pathway.GetRxnsDirs()
        
        curr_max = -1000
        
        pw_r_map = {}
        pw_r_sum = 0
        pw_n_r = 0
        
        if (len(rxns_dict) < 2):
            n_short_pw += 1
        else:
            has_thermo_data = 0
            for rxn,pos in rxns_dict.iteritems():
                if rxn not in total_rxns_map:
                    total_rxns_map[rxn] = 1
                    total_rxns += 1
                    
                try:
                    reg_type = metacyc.GetRegulationType(rxn)
                    if reg_type != None:
                        if reg_type not in reg_hist:
                            reg_hist[reg_type] = [0] * 100
                        reg_hist[reg_type][pos] += 1

                    sparse = metacyc.rxn_uid2sparse_reaction(rxn)
                    
                    if (not sparse):
                        #print 'Reaction not found, uid: %s' % rxn
                        if rxn not in rxns_map:
                            rxns_map[rxn] = 1
                            rxn_map_misses += 1
                        continue
                    
                    sparse = kegg.BalanceReaction(metacyc.sparse2kegg_cids(sparse, kegg), balance_water=True)
                    
                    if (rxn in rxn_dirs):
                        flux = rxn_dirs[rxn]
                    else:
                        if rxn not in rxns_map:
                            rxns_map[rxn] = 1
                            flux_misses += 1
                        continue
                    
                    # deltaG r = CalculateReversabilitydeltaG(sparse, thermo, c_mid, pH, pMg, I, T, concentration_map=cmap)
                    r = CalculateReversability(sparse, thermo, c_mid, pH, pMg, I, T,
                                                      concentration_map=cmap)
                    
                    has_thermo_data = 1
                    
                    if (r == None):
                        if rxn not in rxns_map:
                            rxns_map[rxn] = 1
                            n_only_currency += 1
                    else:
                        r *= flux
                        #r = abs(r)
                                                
                        pw_n_r += 1
                        pw_r_sum += abs(r)
                        if ((pos - 1) in pw_r_map):
                            pw_r_map[pos - 1] = pw_r_map[pos - 1] + ';' + str(r)
                        else:
                            pw_r_map[pos - 1] = str(r)
                            
                        debug_file.write("%s\t%d\t%s\t%s\t%s\t%f\n" % (str(pathway.name), pos, metacyc.uid2reaction_map[rxn].name, metacyc.uid2reaction_map[rxn].ec_number, metacyc.uid2reaction_map[rxn].equation, r))
                        
                        if pos >= 2:
                            histogram['Not first'].append(r)
                        if pos >= 6:
                            histogram['Rest'].append(r)
                        else:
                            histogram.setdefault(pos, []).append(r)
                        if rxn not in rxns_map:
                            rxns_map[rxn] = 1
                            hits += 1
                            
                except thermodynamics.MissingCompoundFormationEnergy:
                    if rxn not in rxns_map:
                        rxns_map[rxn] = 1
                        misses += 1
                    continue
                except KeggNonCompoundException:
                    if rxn not in rxns_map:
                        rxns_map[rxn] = 1
                        comp_translation_misses += 1
                    continue
                except MetaCycNonCompoundException:
                    if rxn not in rxns_map:
                        rxns_map[rxn] = 1
                        meta_compound_misses += 1
                    continue
                except KeggReactionNotBalancedException:
                    if rxn not in rxns_map:
                        rxns_map[rxn] = 1
                        meta_non_balanced += 1
                        #print 'Reaction cannot be balanced, uid: %s' % rxn
                    continue
            
            if has_thermo_data == 1:
                n_pathways_with_thermo_data += 1
                
            if (pw_n_r >= 1):
                avg_abs_r = pw_r_sum / pw_n_r
                if (len(pw_r_map) > 1 and 0 in pw_r_map):
                    first_r = -1000
                    curr_max = -1000
                    n_valid += 1
                    for pos,r_list in pw_r_map.iteritems():
                        for r in r_list.split(';'):
                            r = float(r)
                            if (pos == 0):
                                first_r = r
                            curr_max = r if r > curr_max else curr_max
                             
                            r = float(r) / avg_abs_r
                            rel_histogram.setdefault(pos, []).append(r)
                            if pos > 1:
                                rel_histogram['Not first'].append(r)
                    if (first_r == curr_max):
                        n_first_max += 1
                
    debug_file.close()
    
    debug_file = open('../res/metacyc_' + id + ('_constrained_' if len(cmap) > 0 else '_non_constrained_') + 'rev_stats.txt', 'w')
    debug_file.write("Pathways with positional data: %d\n" % len(metacyc.uid2pathway_map))
    debug_file.write("Pathways shorter then 2: %d\n" % n_short_pw)
    debug_file.write("Pathways with thermodynamic data: %d\n" % n_pathways_with_thermo_data)
    debug_file.write("Total number of reactions: %d\n" % total_rxns)
    debug_file.write("Reactions with known dG0: %d\n" % hits)
    debug_file.write("Reactions with unknown dG0: %d\n" % misses)
    debug_file.write("Reactions with unknown compounds (translation to Kegg failed): %d\n" % comp_translation_misses)
    debug_file.write("Reactions with missing flux direction: %d\n" % flux_misses)
    debug_file.write("Reactions with compounds not found in MetaCyc compounds file: %d\n" % rxn_map_misses)
    debug_file.write("Reactions with compounds not in MetaCyc compounds file: %d\n" % meta_compound_misses)
    debug_file.write("Reactions which couldn't be balanced: %d\n" % meta_non_balanced)
    debug_file.write("Reactions with only currency compounds: %d\n" % n_only_currency)
    debug_file.close()
    
    return histogram,rel_histogram,(float(n_first_max)/n_valid),reg_hist

def calculate_metacyc_regulation_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, metacyc, cmap, id):
    kegg = Kegg.getInstance()
    histogram = {}
    
    
    hits = 0
    misses = 0

    comp_translation_misses = 0
    
    rxn_map_misses = 0
    meta_compound_misses = 0
    meta_non_balanced = 0
    n_only_currency = 0
    
    debug_file = open('../res/metacyc' + id + ('_constrained_' if len(cmap) > 0 else '_non_constrained_') + 'reg.txt', 'w')
    debug_file.write("Regulation type\tReaction Name\tEC list\tEquation\tRev IND\n")

    for rxn_id,rxn in metacyc.uid2reaction_map.iteritems():
        try:
            sparse = metacyc.rxn_uid2sparse_reaction(rxn_id)
                    
            if not sparse:
                #print 'Reaction not found, uid: %s' % rxn
                rxn_map_misses += 1
                continue
                    
            sparse = kegg.BalanceReaction(metacyc.sparse2kegg_cids(sparse, kegg), balance_water=True)            
                    
            r = CalculateReversability(sparse, thermo, c_mid, pH, pMg, I, T,
                                                      concentration_map=cmap)
            if (r == None):
                n_only_currency += 1
            else:
                r = abs(r)
                reg_type = metacyc.GetRegulationType(rxn_id)
                
                debug_file.write("%s\t%s\t%s\t%s\t%f\n" % (reg_type, rxn.name, rxn.ec_number, rxn.equation, r))
                
                histogram.setdefault(reg_type, []).append(r)      
                hits += 1
        except thermodynamics.MissingCompoundFormationEnergy:
            misses += 1
            continue
        except KeggNonCompoundException:
            comp_translation_misses += 1
        except MetaCycNonCompoundException:
            meta_compound_misses += 1
            continue
        except KeggReactionNotBalancedException:
            meta_non_balanced += 1
            #print 'Reaction cannot be balanced, uid: %s' % rxn
            continue
                
    debug_file.close()
    
    debug_file = open('../res/metacyc' + id + ('_constrained_' if len(cmap) > 0 else '_non_constrained_') + 'reg_stats.txt', 'w')
    debug_file.write("Reactions with known dG0: %d\n" % hits)
    debug_file.write("Reactions with unknown dG0: %d\n" % misses)
    debug_file.write("Reactions with unknown compounds (translation to Kegg failed): %d\n" % comp_translation_misses)
    debug_file.write("Reactions with compounds not found in MetaCyc compounds file: %d\n" % rxn_map_misses)
    debug_file.write("Reactions with compounds not in MetaCyc compounds file: %d\n" % meta_compound_misses)
    debug_file.write("Reactions which couldn't be balanced: %d\n" % meta_non_balanced)
    debug_file.write("Reactions with only currency compounds: %d\n" % n_only_currency)
    debug_file.close()
    
    return histogram

    
def analyse_reversibility(thermo, name):
    html_fname = '../res/' + name + '_reversibility.html'
    logging.info('Writing HTML output to %s', html_fname)
    html_writer = HtmlWriter(html_fname)
    c_mid = 1e-4
    cmap = GetConcentrationMap()
    pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)
    
    histogram, rel_histogram, perc_first_max = calculate_reversibility_histogram(
        thermo, c_mid, pH, pMg, I, T, cmap=cmap, id=name)
    
    html_writer.write('<h1>' + name + ': Constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    # deltaG plot fig1 = plot_histogram(histogram, html_writer, title='With constraints on co-factors', legend_loc='lower right' , xlim=80)
    fig1 = plot_histogram(histogram, html_writer, title='With constraints on co-factors', xlim=10)
    
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility1.png', figure=fig1, format='png')
    
    #fig1_bs = plot_bootstrap_stats(histogram, title='With constraints on co-factors')
    #html_writer.embed_matplotlib_figure(fig1_bs, width=640, height=480)
    #pylab.savefig('../res/' + name + '_kegg_reversibility1_bs.png', figure=fig1_bs, format='png')
    
    fig1_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module with constraints on co-factors', xlim=5)
    html_writer.embed_matplotlib_figure(fig1_rel, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility1_rel.png', figure=fig1_rel, format='png')
    
    histogram, rel_histogram, perc_first_max = calculate_reversibility_histogram(
        thermo, c_mid, pH, pMg, I, T, cmap={}, id=name)

    html_writer.write('<h1>' + name + ': Non constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig2 = plot_histogram(histogram, html_writer, title='No constraints on co-factors', xlim=20)
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2.png', figure=fig2, format='png')
    
    fig2_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module, no constraints on co-factors', xlim=5)
    html_writer.embed_matplotlib_figure(fig2_rel, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2_rel.png', figure=fig2_rel, format='png')

    #(histogram,rel_histogram,perc_first_max) = calculate_reversibility_histogram(thermo, c_mid, pH, pMg, I, T,
    #                                              cmap=GetFullConcentrationMap(thermo), id=(name + '_full'))
    
    #html_writer.write('<h1>' + name + ': Constrained metabolites</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    #fig3 = plot_histogram(histogram, html_writer, title='With constraints all metabolites with known concentration', xlim=20)
    #html_writer.embed_matplotlib_figure(fig3, width=640, height=480)
    #pylab.savefig('../res/' + name + '_kegg_reversibility3.png', figure=fig3, format='png')
    
    #fig3_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module with constrainted metabolites', xlim=5)
    #html_writer.embed_matplotlib_figure(fig3_rel, width=640, height=480)
    #pylab.savefig('../res/' + name + '_kegg_reversibility3_rel.png', figure=fig3_rel, format='png')
    
def metacyc_data(org, id, thermo, max_pathway_length_for_fig=8):
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/' + org + '_' + id + '_reversibility.html')
    kegg = Kegg.getInstance()
    metacyc_inst = MetaCyc(org, db)
    c_mid = 1e-4
    cmap = GetConcentrationMap()
    pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)
    
    (histogram,rel_histogram,perc_first_max, reg_hist) = calculate_metacyc_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, metacyc_inst,
                                                  cmap=cmap, id=(org + '_' + id))
    
    html_writer.write('<h1>Constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    # deltaG plot fig1 = plot_histogram(histogram, html_writer, title=('%s pathways: With constraints on co-factors' % org), legend_loc='lower right' , xlim=80)
    fig1 = plot_histogram(histogram, html_writer, title=('%s pathways: With constraints on co-factors' % org), xlim=10)
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id +  '_reversibility1.png', figure=fig1, format='png')

    fig1_reg = plot_bars(reg_hist, title=('%s pathways: Position of regulated reactions' % org), max_pathway_length=max_pathway_length_for_fig)
    html_writer.embed_matplotlib_figure(fig1_reg, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id +  '_reg_rxns.png', figure=fig1_reg, format='png')
    
    fig1_reg_stacked = plot_stacked_bars(reg_hist, title=('%s pathways: Position of regulated reactions' % org), max_pathway_length=max_pathway_length_for_fig)
    html_writer.embed_matplotlib_figure(fig1_reg_stacked, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id +  '_reg_stacked_rxns.png', figure=fig1_reg_stacked, format='png')

    #fig1_bs = plot_bootstrap_stats(histogram, title=('%s pathways: With constraints on co-factors' % org))
    #html_writer.embed_matplotlib_figure(fig1_bs, width=640, height=480)
    #pylab.savefig('../res/' + org + '_' + id +  '_reversibility1_bs.png', figure=fig1_bs, format='png')    
    

    fig1_rel = plot_histogram(rel_histogram, html_writer, title=('%s pathways: Normed per pathway with constraints on co-factors' % org), xlim=5)
    html_writer.embed_matplotlib_figure(fig1_rel, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id + '_reversibility1_rel.png', figure=fig1_rel, format='png')
    
    (histogram,rel_histogram,perc_first_max, reg_hist) = calculate_metacyc_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, metacyc_inst,
                                                  cmap={}, id=(org + '_' + id))
    
    html_writer.write('<h1>Non constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig2 = plot_histogram(histogram, html_writer, title=('%s pathways: No constraints on co-factors' % org ), xlim=20)
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id + '_reversibility2.png', figure=fig1, format='png')

    fig2_rel = plot_histogram(rel_histogram, html_writer, title=('%s pathways: Normed per pathway no constraints on co-factors' % org), xlim=5)
    html_writer.embed_matplotlib_figure(fig2_rel, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + id + '_reversibility2_rel.png', figure=fig2_rel, format='png')

def meta_regulated_rxns_cumul_plots(org, id, thermo):
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/' + org + id + '_regulation.html')
    kegg = Kegg.getInstance()
    metacyc_inst = MetaCyc(org, db)
    c_mid = 1e-4
    cmap = GetConcentrationMap()
    pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)
    
    histogram = calculate_metacyc_regulation_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, metacyc_inst,
                                                  cmap=cmap, id=(org + id))
    
    html_writer.write('<h1>Constrained co-factors</h1>')
    fig1 = plot_histogram(histogram, html_writer, title=('%s Reactions: With constraints on co-factors' % org), xlim=20, min_to_show=5, xmin=0, legend_loc='lower right')
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + org + id +  '_regulation.png', figure=fig1, format='png')

def calc_palsson_rxns_irrev(org, fig_name, rxns_file, comp2cid_file, thermo, c_mid, pH, pMg, I, T, cmap, xlim=40):
    html_fname = '../res/' + org + '_' + fig_name + '_model_rev.html'
    logging.info('Writing HTML output to %s', html_fname)
    html_writer = HtmlWriter(html_fname)
    kegg = Kegg.getInstance()

    histogram = {}
    
    comp2cid_map = {}
    map_file = open (comp2cid_file, 'r')
    line = map_file.readline().strip()
    
    n_miss_comp = 0
    n_not_balanced = 0
    n_hits = 0
    misses = 0
    n_only_currency = 0
    
    debug_file = open('../res/fba_' + org + '_' + fig_name + '_rev.txt', 'w')
    debug_file.write("Reaction Name\tFBA classification\tEquation\tRev IND\n")

    # Read compounds ids 2 Kegg cid's mapping file into a dict
    while (line):
        if (len(line) == 0):
            continue
        (id,cid) = line.split('\t')
        comp2cid_map[id] = cid
        line = map_file.readline().strip()
    map_file.close()

    left = ''
    right = ''

    rxns_input = open(rxns_file, 'r')
    line = rxns_input.readline().strip()
    
    while (line):
        line = rxns_input.readline().strip() # Skipping the header line
        if (len(line) == 0):
            continue
        
        fields = line.split('\t')
        name = fields[1]
        equation = fields[3]
        reversible = fields[7]
        
        # All fields: abbrev, name, syn, equation, subsystem, compartment, ecnumber, reversible, translocation, internal_id, confidence_score, notes 

        sparse = {}
        equation = re.sub('\[[pce]\]', '', equation)
        equation = re.sub(' : ', '', equation)
        left,right = re.split('<==>|-->', equation)
        
        try:
            parse_palsson_side_eq (left, -1, sparse, comp2cid_map)
            parse_palsson_side_eq (right, 1, sparse, comp2cid_map)
        
            sparse = kegg.BalanceReaction(sparse, balance_water=True)
        
            r = CalculateReversability(sparse, thermo, c_mid, pH, pMg, I, T,
                                                      concentration_map=cmap)
            if r == None:
                n_only_currency +=1
                continue
            else:
                #r = abs(r)
                n_hits += 1
                histogram.setdefault(reversible, []).append(r)
                debug_file.write('%s\t%s\t%s\t%f\n' % (name, reversible, equation, r))
        except KeggNonCompoundException:
            logging.debug('No Kegg cid for reaction %s' % name)
            n_miss_comp += 1
            continue
        except KeggReactionNotBalancedException:
            logging.debug('Reaction %s is not balanced' % name)
            n_not_balanced += 1
            continue
        except thermodynamics.MissingCompoundFormationEnergy:
            misses += 1
            continue
        except IrrevParseException:
            logging.debug('Error parsing reaction %s' % name)
            continue
        
    debug_file.close() 
    
    debug_file = open('../res/fba_' + org + '_' + fig_name + '_rev_stats.txt', 'w')
    debug_file.write("Reactions with known dG0: %d\n" % n_hits)
    debug_file.write("Reactions with unknown dG0: %d\n" % misses)
    debug_file.write("Reactions with unknown compounds (translation to Kegg failed): %d\n" % n_miss_comp)
    debug_file.write("Non balanced reactions: %d\n" % n_not_balanced)
    debug_file.write("Reactions of solely currency copmounds: %d\n" % n_only_currency)
    debug_file.close()
    
    # plot the bar 
    fig = pylab.figure()
    pylab.hold(True)

    colors = {'Reversible':'blue', 'Irreversible':'red'}
    
    i = 1
    for key, value in histogram.iteritems():
        new_value = map(lambda x: min(x, xlim), value)
        new_value = map(lambda x: max(x, xlim * - 1), new_value)
                 
        pylab.subplot (2,1,i)
        i +=1
        label = '%s (%d reactions)' % (key, len(new_value))      
        #pylab.hist(value, 101, range=(-50,50), normed=True, label=label, color=colors[key])
        pylab.hist(new_value, xlim * 2 + 1, range=(xlim * -1 ,xlim), normed=True, label=label, color=colors[key])
        legendfont = matplotlib.font_manager.FontProperties(size=9)
        pylab.legend(prop=legendfont)
        pylab.xlabel('irreversability')
        pylab.ylabel('Fraction of reactions')
        #pylab.title(org + '_' + id + ' model')
    
    pylab.hold(False)
    
    html_writer.write('<h1>Irreversibility index histogram on Palsson\'s model %s</h1>' % org)
    html_writer.embed_matplotlib_figure(fig, width=640, height=480)
    pylab.savefig('../res/' + org + '_' + fig_name + '_model_irrev.png', figure=fig, format='png')
    
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

class IrrevParseException(Exception):
    pass

def calc_cons_rxns_corr(thermo, name):
    html_fname = '../res/' + name + '_rev_pair_corr.html'
    logging.info('Writing HTML output to %s', html_fname)
    html_writer = HtmlWriter(html_fname)
    c_mid = 1e-4
    cmap = GetConcentrationMap()
    pH, pMg, I, T = (7.0, 3.0, 0.25, 298.15)
    
    (first, second) = get_reversibility_consecutive_pairs(thermo, c_mid, pH, pMg, I, T,
                                                  cmap=cmap, id=name)
    html_writer.write('<h1>' + name + ': Constrained co-factors</h1><br>')
    
    fig1 = cons_pairs_dot_plot (first, second, xlim=50)
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + name + '_rev_pairs_corr.png', figure=fig1, format='png')
    
    fig2 = cons_pairs_dot_plot (first, second, xlim=10)    
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + name + '_rev_pairs_corr_zoom.png', figure=fig2, format='png')

def cons_pairs_dot_plot(first, second, xlim=30):
    fig = pylab.figure()
    pylab.hold(True)
        
    pylab.plot(first, second, marker='.', linestyle='None', markersize=5)
    
    r2 = util.calc_r2(first,second)
    pearson = numpy.corrcoef(first, second)[0,1]
    pylab.title('Consecutive reactions irreversibility (R^2=%.2f, Pearson=%.2f)' % (r2, pearson), fontsize=14)
    
    pylab.xlabel('First reaction: irreversibility', fontsize=12)
    pylab.ylabel('Second reaction: irreversibility', fontsize=12)
    min_x = min(first)
    max_x = max(first)
    pylab.plot([min_x, max_x], [min_x, max_x], 'r--')
    pylab.plot([min_x, max_x], [0,0.001], 'k-')
    pylab.plot([0,0.001], [min_x, max_x], 'k-')
    min_x = max(min_x, -1 * xlim)
    max_x = min(max_x, xlim)
    pylab.axis([min_x, max_x, min_x, max_x])
    
    pylab.hold(False)
    
    return fig
 
def get_reversibility_consecutive_pairs(G, c_mid, pH, pMg, I, T, cmap, id):
    first = []
    second = []

    hits = 0
    misses = 0
    only_currency = 0
    non_balanced = 0
    kegg = Kegg.getInstance()

    fname = '../res/kegg_' + id + '_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev_pairs.txt'
    logging.info('Writing output to %s', fname)
    debug_file = open(fname, 'w')
    debug_file.write('Module\tRxn 1\tIrrev 1\tRxn 2\tIrrev 2\n')
        
    for mid, rid_flux_list in kegg.mid2rid_map.iteritems():
        if not rid_flux_list or len(rid_flux_list) < 2:
            continue
        
        prev_r = None
        prev_i = -1
        prev_name = None
        

        
        for i, (rid, flux) in enumerate(rid_flux_list):
            try:
                sparse = kegg.rid2sparse_reaction(rid)
                Reaction.BalanceSparseReaction(sparse, balance_water=True)
                r = CalculateReversability(sparse, G, c_mid, pH, pMg, I, T,
                                           concentration_map=cmap)
                if r == None:
                    only_currency += 1
                    continue
                
                r *= flux
                         
                if prev_r != None and i - prev_i == 1:
                    first.append(prev_r)
                    second.append(r)
                    debug_file.write('%s\t%s\t%f\t%s\t%f\n' % (re.sub('\t', ' ', kegg.mid2name_map[mid]), prev_name, prev_r, kegg.rid2name(rid), r ))
                
                prev_r = r
                prev_i = i
                prev_name = kegg.rid2name(rid)
                hits += 1                
                            
            except thermodynamics.MissingCompoundFormationEnergy:
                misses += 1
                continue
            except KeggReactionNotBalancedException:
                    non_balanced += 1
                    #print 'Reaction cannot be balanced, uid: %s' % rxn
                    continue
            
    
    debug_file.close()
    
    debug_file = open('../res/kegg_' + id + '_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev_pairs_corr_stats.txt', 'w')
    debug_file.write("Reactions with known dG0: %d\n" % hits)
    debug_file.write("Reactions with unknown dG0: %d\n" % misses)
    debug_file.write("Non balanced reactions: %d\n" % non_balanced)
    debug_file.write("Reactions without unknown concentrations: %d\n" % only_currency)
    debug_file.close()

    return first,second

def compare_reversibility_to_dG0(thermo, name):
    html_fname = '../res/' + name + '_rev_vs_dG.html'
    logging.info('Writing HTML output to %s', html_fname)
    html_writer = HtmlWriter(html_fname)
    kegg = Kegg.getInstance()
    c_mid = 1e-4
    pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)
    
    x_range = (1e-9, 1e9)
    y_range = (1e-9, 1e9)

    x_threshold = 1e3; x_color = 'blue'
    y_threshold = 1e3; y_color = 'red'
    
    # plot the profile graph
    pylab.rcParams['text.usetex'] = False
    pylab.rcParams['legend.fontsize'] = 10
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 14
    pylab.rcParams['lines.linewidth'] = 2
    pylab.rcParams['lines.markersize'] = 6
    pylab.rcParams['figure.figsize'] = [6.0, 6.0]
    pylab.rcParams['figure.dpi'] = 100
    
    fig = pylab.figure()
    pylab.xlabel(r"$K^'$", figure=fig)
    pylab.ylabel(r"$\gamma = \left( K^' / \Gamma(1) \right)^{2/N}$", figure=fig)
    pylab.axvspan(x_range[0], 1.0/x_threshold, ymin=0, ymax=1, color=x_color, alpha=0.3)
    pylab.axvspan(x_threshold, x_range[1], ymin=0, ymax=1, color=x_color, alpha=0.3)
    pylab.axhspan(y_range[0], 1.0/y_threshold, xmin=0, xmax=1, color=y_color, alpha=0.3)
    pylab.axhspan(y_threshold, y_range[1], xmin=0, xmax=1, color=y_color, alpha=0.3)

    stoichiometries = [(1, 1, 'orange'), 
                       (1 ,2, 'r--'), 
                       (2, 1, 'g--'), 
                       (2, 2, 'blue'), 
                       (2, 3, 'm--'), 
                       (3, 2, 'c--'), 
                       (3, 3, 'pink')]
    fig.hold(True)
    for n_s, n_p, style in stoichiometries:
        gamma = [(Keq*c_mid**(n_p - n_s)) ** (2.0/(n_p + n_s)) for Keq in x_range]
        pylab.plot(x_range, gamma, style, figure=fig, label="%d:%d" % (n_s, n_p))
    pylab.legend(loc='upper left')
    #reactions = []
    #reactions.append(kegg.rid2reaction[1068]) # aldolase
    #reactions.append(kegg.rid2reaction[24]) # aldolase
    #for reaction in
    reactions = kegg.AllReactions() 
    
    counters = {}
    data_mat = pylab.zeros((0, 4))
    for reaction in reactions:
        #if reaction.rid % 10 != 1:
        #    continue
        try:
            reaction.Balance(balance_water=True)

            dG0 = reaction.PredictReactionEnergy(thermo, pH, pMg, I, T) 
            Keq = pylab.exp(-dG0/(R*T))

            n_s = -sum([x for cid, x in reaction.sparse.iteritems() if (x < 0 and cid != 1)])
            n_p = sum([x for cid, x in reaction.sparse.iteritems() if (x > 0 and cid != 1)])
            if (n_p + n_s) == 0:
                continue
            
            gamma = (Keq*c_mid**(n_p - n_s)) ** (2.0/(n_p + n_s))
            
            if Keq < 1.0/x_threshold:
                Krev = -1
            elif Keq < x_threshold:
                Krev = 0
            else:
                Krev = 1
                
            if gamma < 1.0/y_threshold:
                Grev = -1
            elif gamma < y_threshold:
                Grev = 0
            else:
                Grev = 1
                
            counters.setdefault((Krev, Grev), 0)
            counters[Krev, Grev] += 1
            data_mat = pylab.vstack([data_mat, [Keq, gamma, Krev, Grev]])
        except (MissingCompoundFormationEnergy, KeggReactionNotBalancedException):
            pass
    
    fig.hold(True)
    for Krev, Grev in counters.keys():
        x_pos = x_threshold ** (Krev*2)
        y_pos = y_threshold ** (Grev*2)
        pylab.text(x_pos, y_pos, "%.1f%%" % (100.0 * counters[Krev, Grev] / data_mat.shape[0]), 
                   horizontalalignment='center',
                   verticalalignment='center')

    pylab.xscale('log', figure=fig)
    pylab.yscale('log', figure=fig)
    pylab.ylim(y_range, figure=fig)
    pylab.xlim(x_range, figure=fig)
    html_writer.embed_matplotlib_figure(fig, width=500, height=500, name=name + "_k_vs_g")
   
    fig = pylab.figure()
    max_gamma = 1e20
    pylab.axvspan(y_threshold, max_gamma, ymin=0, ymax=1, color=y_color, alpha=0.3)
    abs_gamma = pylab.exp(abs(pylab.log(data_mat[:,1])))
    plotting.cdf(abs_gamma, label='gamma', figure=fig)
    pylab.xscale('log')
    pylab.xlim((1, max_gamma))
    html_writer.embed_matplotlib_figure(fig, width=500, height=500, name=name + "_cdf")

def main():
    #logging.getLogger('').setLevel(logging.DEBUG)
    #try_kegg_api()
    estimators = LoadAllEstimators()
    #analyse_reversibility(estimators['hatzi_gc'], 'HatziGC')
    #analyse_reversibility(estimators['milo_gc'], 'MiloGC_zoom')
    
    compare_reversibility_to_dG0(estimators['milo_gc'], 'MiloGC')
    #compare_reversibility_to_dG0(estimators['alberty'], 'Alberty')
    #metacyc_data('meta', 'zoom', estimators['milo_gc'], max_pathway_length_for_fig=6)
    #metacyc_data('ecoli', 'zoom', estimators['milo_gc'], max_pathway_length_for_fig=4)
    
    #meta_regulated_rxns_cumul_plots('meta', '_abs', estimators['milo_gc'])
    #meta_regulated_rxns_cumul_plots('ecoli', '_abs', estimators['milo_gc'])
    
    #calc_palsson_rxns_irrev('ecoli', 'zoom', '../data/reactionList_iAF1260.txt',
    #                        '../data/ecoli2kegg_cid.txt',
    #                        thermo=estimators['milo_gc'], c_mid=1e-3,
    #                        pH=7.0, pMg=3.0, I=0.1, T=298.15,
    #                        cmap=GetConcentrationMap(Kegg.getInstance()),
    #                        xlim=20)
    #calc_cons_rxns_corr(estimators['milo_gc'], 'MiloGC')
    
if __name__ == "__main__":
    main()