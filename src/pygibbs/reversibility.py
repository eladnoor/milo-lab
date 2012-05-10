#!/usr/bin/python

from scipy import stats
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs import thermodynamics, thermodynamic_comparison,\
    thermodynamic_constants
from pygibbs.thermodynamic_constants import R
from pygibbs.groups import GroupContribution
import pylab
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggNonCompoundException
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.kegg_reaction import Reaction
from toolbox import plotting
from SOAPpy import WSDL 
from pygibbs.metacyc import MetaCyc, MetaCycNonCompoundException
import logging
import matplotlib
import re
import numpy as np
import random
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy
from pygibbs.thermodynamic_errors import MissingReactionEnergy
from pygibbs.feist_ecoli import Feist
from pygibbs.compound_abundance import CompoundAbundance
from toolbox.molecule import OpenBabelError

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
            reaction = G.kegg.rid2reaction(rid)
            r = CalculateReversability(reaction, G, c_mid, pH, pMg, I, T)
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
DEFAULT_CMID = 1e-4
DEFAULT_PH, DEFAULT_PMG, DEFAULT_I, DEFAULT_T = (7.0, 3.0, 0.1, 298.15)

def GetEmptyConcentrationMap():
    return {WATER:1, HPLUS:1}

def GetConcentrationMap():
    kegg = Kegg.getInstance()
    cmap = GetEmptyConcentrationMap() 
    for cid in kegg.get_all_cids():
        lower, upper = kegg.get_bounds(cid)
        if lower and upper:
            # In the file we got this data from lower = upper 
            cmap[cid] = lower
    return cmap

def GetFullConcentrationMap():
    abundance = CompoundAbundance.LoadConcentrationsFromBennett()
    cmap = GetEmptyConcentrationMap() 
    kegg = Kegg.getInstance()
    
    for cid in kegg.get_all_cids():
        cmap[cid] = abundance.GetConcentration(cid, c0=None, medium='glucose')
        if cmap[cid] == None:
            lower, upper = kegg.get_bounds(cid)
            if lower and upper:
                # In the file we got this data from lower = upper 
                cmap[cid] = lower
            else:
                cmap.pop(cid)
    return cmap

def ConcentrationFactor(reaction,
                        concentration_map,
                        c_mid):
    factor = 0.0
    for cid, stoic in reaction.sparse.iteritems():
        concentration = concentration_map.get(cid, None) or c_mid
        factor += pylab.log(concentration) * stoic
    return factor

def CalculateReversability(reaction, thermo, concentration_map=None, logscale=False):
    cmap = concentration_map or GetEmptyConcentrationMap()
    dG0 = reaction.PredictReactionEnergy(thermo)
    ln_Gamma = ConcentrationFactor(reaction, cmap, thermo.c_mid)
    sum_abs_s = sum([abs(x) for k, x in reaction.sparse.iteritems()
                     if k not in cmap])
    if sum_abs_s == 0:
        return None
    else:
        if logscale: # 2/N * (log K' - log Q'')
            return 2.0 / sum_abs_s * (-dG0/(R*thermo.T) - ln_Gamma) 
        else:        # (K' / Q'') ^ (2/N)
            return np.exp(-dG0/(R*thermo.T) - ln_Gamma) ** (2.0 / sum_abs_s)

def CalculateReversabilityV2(reaction, thermo, concentration_map=None):
    cmap = concentration_map or GetEmptyConcentrationMap()
    dG0 = reaction.PredictReactionEnergy(thermo)
    cfactor = ConcentrationFactor(reaction, cmap, thermo.c_mid)
    
    if (cfactor == 0):
        return None
    else:
        return pylab.exp(-dG0/(R * thermo.T * cfactor))

def CalculateReversabilitydeltaG(reaction, thermo, concentration_map=None):
    cmap = concentration_map or GetEmptyConcentrationMap()
    dG0 = reaction.PredictReactionEnergy(thermo)
    cfactor = ConcentrationFactor(reaction, cmap, thermo.c_mid)
    return pylab.exp(dG0 + R * thermo.T * cfactor)

def calculate_reversibility_histogram(G, thermo, cmap, id):
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
                gamma = CalculateReversability(reaction, thermo, 
                                               concentration_map=cmap)
                r = pylab.log10(gamma)
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
                dbg = "%s__DEL__%d__DEL__%s__DEL__%s__DEL__%s__DEL__%s__DEL__%f\n" % (kegg.mid2name_map[mid], i+1, rxn.name, rxn.definition, str(rxn.ec_list), rxn.equation, r)
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

def plot_bars(pos_count, title='', max_pathway_length=8, legend_loc='upper right'):
    n_labels = len(pos_count)
    ind = np.arange(max_pathway_length)
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

def plot_stacked_bars(pos_count, title='', max_pathway_length=8, legend_loc='upper right'):
    
    ind = np.arange(max_pathway_length)
    width = 0.35
    
    fig = pylab.figure()
    pylab.hold(True)
    ax = fig.add_subplot(111)

    colors = {'No known regulation':'grey', 'Activated':'green', 'Inhibited':'red', 'Mixed regulation':'blue'}
    plot_order = ['Inhibited', 'Mixed regulation', 'Activated', 'No known regulation']
    
    total = np.array([0] * max_pathway_length)
    for label,vals in pos_count.iteritems():
        total += np.array(vals[1:max_pathway_length+1])

    so_far = np.array([0] * max_pathway_length)
    for label in plot_order:
        vals = pos_count[label]
        curr_vals = np.array(vals[1:max_pathway_length+1]) * 1.0 / total
        ax.bar(ind, tuple(curr_vals), width, color=colors[label], bottom=tuple(so_far), label=('%s (%d)' % (label, sum(vals[1:max_pathway_length+1]))))
        so_far = so_far + np.array(curr_vals)
    
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



def calculate_metacyc_reversibility_histogram(thermo, metacyc, cmap, id):
    
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
            for rxn, pos in rxns_dict.iteritems():
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
                    
                    reaction = Reaction(str(rxn), sparse=metacyc.sparse2kegg_cids(sparse, kegg))
                    
                    if (rxn in rxn_dirs):
                        flux = rxn_dirs[rxn]
                    else:
                        if rxn not in rxns_map:
                            rxns_map[rxn] = 1
                            flux_misses += 1
                        continue
                    
                    # deltaG r = CalculateReversabilitydeltaG(sparse, thermo, c_mid, pH, pMg, I, T, concentration_map=cmap)
                    gamma = CalculateReversability(reaction, thermo,
                                                   concentration_map=cmap)
                    r = pylab.log10(gamma)
                    
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

def calculate_metacyc_regulation_reversibility_histogram(thermo, metacyc, cmap, id):
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
            reaction = Reaction(rxn_id, metacyc.sparse2kegg_cids(sparse, kegg))
            reaction.Balance(balance_water=True)            
            gamma = CalculateReversability(reaction, thermo,
                                           concentration_map=cmap)
            if gamma == None:
                n_only_currency += 1
            else:
                gamma = max(gamma, 1/gamma)
                reg_type = metacyc.GetRegulationType(rxn_id)
                
                debug_file.write("%s\t%s\t%s\t%s\t%f\n" % (reg_type, rxn.name, rxn.ec_number, rxn.equation, gamma))
                
                histogram.setdefault(reg_type, []).append(gamma)      
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
    cmap = GetConcentrationMap()
    
    histogram, rel_histogram, perc_first_max = calculate_reversibility_histogram(
        thermo, cmap=cmap, id=name)
    
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
        thermo, cmap={}, id=name)

    html_writer.write('<h1>' + name + ': Non constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig2 = plot_histogram(histogram, html_writer, title='No constraints on co-factors', xlim=20)
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2.png', figure=fig2, format='png')
    
    fig2_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module, no constraints on co-factors', xlim=5)
    html_writer.embed_matplotlib_figure(fig2_rel, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2_rel.png', figure=fig2_rel, format='png')

    #(histogram,rel_histogram,perc_first_max) = calculate_reversibility_histogram(thermo, c_mid, pH, pMg, I, T,
    #                                              cmap=GetFullConcentrationMap(), id=(name + '_full'))
    
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
    metacyc_inst = MetaCyc(org, db)
    cmap = GetConcentrationMap()
    
    (histogram,rel_histogram,perc_first_max, reg_hist) = calculate_metacyc_reversibility_histogram(thermo, metacyc_inst,
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
    
    (histogram,rel_histogram,perc_first_max, reg_hist) = calculate_metacyc_reversibility_histogram(thermo, metacyc_inst,
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

def compare_annotations(reaction_list, thermo, html_writer, cmap, xlim=1e9):
    html_writer.write('<h1>Compare reaction annotations to Reversibility Index</h1>\n')
    histogram = {}
    error_counts = {'hits': 0, 'misses': 0, 'no_gamma': 0}

    debug_dict_list = []
    for reaction in reaction_list:
        try:
            dG0 = reaction.PredictReactionEnergy(thermo)
        except (MissingCompoundFormationEnergy, MissingReactionEnergy) as e:
            logging.warning(str(e))
            error_counts['misses'] += 1
            continue

        gamma = CalculateReversability(reaction, thermo, concentration_map=cmap)
        if gamma is None:
            error_counts['no_gamma'] += 1
        else:
            error_counts['hits'] += 1
            histogram.setdefault(reaction.direction, []).append(gamma)
            debug_dict_list.append({'sortkey':gamma,
                                    'Reaction Name':reaction.name,
                                    'annotation':reaction.direction,
                                    'KEGG Reaction':reaction.to_hypertext(),
                                    'Rev. index':"%.3g" % gamma,
                                    'dG0':"%.2f" % dG0})
    
    debug_dict_list.sort(key=lambda(x):x['sortkey'])
    div_id = html_writer.insert_toggle()
    html_writer.div_start(div_id)
    html_writer.write_table(debug_dict_list, headers=['Rev. index', 
        'dG0', 'Reaction Name', 'KEGG Reaction', 'annotation'])
    html_writer.div_end()
    html_writer.write('</br>\n')
    
    html_writer.write_ul(["Reactions with known dG0: %d" % error_counts['hits'],
                          "Reactions with unknown dG0: %d" % error_counts['misses'],
                          "Reactions with unknown gamma: %d" % error_counts['no_gamma']])
    
    # plot the bar 
    fig = pylab.figure(figsize=(6,6), dpi=90)
    pylab.hold(True)
    plotting.cdf(histogram['<=>'], label='reversible (%d reactions)' % len(histogram['<=>']),
                 style='green', figure=fig)
    plotting.cdf(histogram['=>'], label='forward only (%d reactions)' % len(histogram['=>']),
                 style='red', figure=fig)
    plotting.cdf(histogram['<='], label='reverse only (%d reactions)' % len(histogram['<=']),
                 style='orange', figure=fig)

    pylab.xlabel('Reversability index - $\hat{\gamma}$')
    pylab.ylabel('Cumulative Distribution')
    pylab.xscale('log')
    pylab.xlim((1/xlim, xlim))
    pylab.legend(loc='upper left')
    html_writer.embed_matplotlib_figure(fig, width=640, height=480, name='FEIST_CDF')
    
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
    
    pearson = np.corrcoef(first, second)[0,1]
    pylab.title('Consecutive reactions irreversibility (Pearson=%.2f)' % (pearson), fontsize=14)
    
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
                reaction = kegg.rid2reaction(rid)
                reaction.Balance(balance_water=True)
                r = CalculateReversability(reaction, G, c_mid, pH, pMg, I, T,
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

def compare_reversibility_to_dG0(reaction_list, thermo, html_writer, cmap=None):
    html_writer.write('<h1>Reversibility index vs. equilibrium constants</h1>\n')
    cmap = cmap or GetEmptyConcentrationMap()
    
    x_range = (1e-9, 1e9)
    y_range = (1e-9, 1e9)

    x_threshold = 1e3
    y_threshold = 1e3
    
    regime_counters = {}
    stoich_counters = {}
    data_mat = pylab.zeros((0, 4))
    
    debug_dict_list = []
    for reaction in reaction_list:
        debug_dict = {'name':reaction.name, 
                      'KEGG Reaction':reaction.to_hypertext()}
        
        try:
            reaction.Balance(balance_water=True, exception_if_unknown=True)
        except (KeggReactionNotBalancedException, OpenBabelError):
            continue
        
        dG0 = reaction.PredictReactionEnergy(thermo)
        if np.isnan(dG0):
            debug_dict['sortkey'] = 0
            debug_dict['error'] = "Cannot calculate Gibbs energy"
        else:
            Keq = pylab.exp(-dG0/(R*thermo.T))
    
            n_s = -sum([x for cid, x in reaction.sparse.iteritems() if (x < 0 and cid not in cmap)])
            n_p = sum([x for cid, x in reaction.sparse.iteritems() if (x > 0 and cid not in cmap)])
            if (n_p + n_s) == 0:
                continue
            stoich_counters.setdefault((n_s, n_p), 0)
            stoich_counters[n_s, n_p] += 1
            
            log_gamma = CalculateReversability(reaction, thermo,
                                               concentration_map=cmap,
                                               logscale=True)
            
            if Keq < 1.0/x_threshold:
                Krev = -1
            elif Keq < x_threshold:
                Krev = 0
            else:
                Krev = 1
                
            if log_gamma < -np.log(y_threshold):
                Grev = -1
            elif log_gamma < np.log(y_threshold):
                Grev = 0
            else:
                Grev = 1
                
            regime_counters.setdefault((Krev, Grev), 0)
            regime_counters[Krev, Grev] += 1
            data_mat = pylab.vstack([data_mat, [Keq, log_gamma, Krev, Grev]])
            debug_dict['sortkey'] = log_gamma
            debug_dict['log(&gamma;)'] = "%.2e" % log_gamma
            debug_dict[thermodynamic_constants.symbol_dr_G0_prime] = dG0
        
        debug_dict_list.append(debug_dict)
    
    debug_dict_list.sort(key=lambda(x):x['sortkey'])
    div_id = html_writer.insert_toggle()
    html_writer.div_start(div_id)
    html_writer.write_table(debug_dict_list, headers=['log(&gamma;)',
        thermodynamic_constants.symbol_dr_G0_prime, 'name', 'KEGG Reaction',
        'error'])
    html_writer.div_end()
    html_writer.write('</br>\n')
    
    fig = pylab.figure(figsize=(6,6), dpi=90)
    pylab.xlabel("$K'$", figure=fig)
    pylab.ylabel(r"$\hat{\gamma} = \left( K' / Q'' \right)^{2/N}$", figure=fig)
    
    shading_color = (1.0, 0.7, 0.7)
    #pylab.axvspan(x_range[0], 1.0/x_threshold, ymin=0, ymax=1, color=x_color, alpha=0.3)
    #pylab.axvspan(x_threshold, x_range[1], ymin=0, ymax=1, color=x_color, alpha=0.3)
    #pylab.axhspan(y_range[0], 1.0/y_threshold, xmin=0, xmax=1, color=y_color, alpha=0.3)
    #pylab.axhspan(y_threshold, y_range[1], xmin=0, xmax=1, color=y_color, alpha=0.3)
    pylab.axvspan(x_range[0], 1.0/x_threshold, ymin=1.0/3.0, ymax=2.0/3.0, color=shading_color)
    pylab.axvspan(x_threshold, x_range[1], ymin=1.0/3.0, ymax=2.0/3.0, color=shading_color)
    pylab.axhspan(y_range[0], 1.0/y_threshold, xmin=1.0/3.0, xmax=2.0/3.0, color=shading_color)
    pylab.axhspan(y_threshold, y_range[1], xmin=1.0/3.0, xmax=2.0/3.0, color=shading_color)

    # draw the lines for the specific reaction stoichiometries
    stoichiometries = [(1, 1, '-', '#e49b1c'), 
                       (1 ,2, '--', '#1ce463'), 
                       (2, 1, '--', '#1d1de3'), 
                       (2, 2, '-', '#e41c63')] 
    fig.hold(True)
    for n_s, n_p, style, color in stoichiometries:
        percent = 100.0 * stoich_counters.get((n_s, n_p), 0) / sum(stoich_counters.values())
        gamma = [(Keq / thermo.c_mid**(n_p - n_s)) ** (2.0/(n_p + n_s)) for Keq in x_range]
        pylab.plot(x_range, gamma, style, color=color, linewidth=3,
                   figure=fig, label="%d:%d (%d%%)" % (n_s, n_p, np.round(percent)))
    pylab.legend(loc='upper left')

    for Krev, Grev in regime_counters.keys():
        x_pos = x_threshold ** (Krev*2)
        y_pos = y_threshold ** (Grev*2)
        pylab.text(x_pos, y_pos, "%.1f%%" % (100.0 * regime_counters[Krev, Grev] / data_mat.shape[0]), 
                   horizontalalignment='center',
                   verticalalignment='center')

    pylab.xscale('log', figure=fig)
    pylab.yscale('log', figure=fig)
    pylab.ylim(y_range)
    pylab.xlim(x_range)
    pylab.xticks([1e-9, 1e-6, 1e-3, 1, 1e3, 1e6, 1e9])
    pylab.yticks([1e-9, 1e-6, 1e-3, 1, 1e3, 1e6, 1e9])
    html_writer.embed_matplotlib_figure(fig, width=400, height=400, name="reversibility_vs_keq")
   
    fig = pylab.figure(figsize=(2,2), dpi=90)
    abs_gamma = np.exp(abs(data_mat[:,1]))
    plotting.cdf(abs_gamma, label='gamma', figure=fig)
    pylab.plot([x_threshold, x_threshold], [0, 1], 'k--', figure=fig)
    pylab.xscale('log', figure=fig)
    #pylab.xlabel(r'$\hat{\gamma}$', figure=fig)
    #pylab.ylabel(r'CDF($\hat{\gamma}$)', figure=fig)
    pylab.text(1e6, 0.4, r'CDF($\hat{\gamma}$)', horizontalalignment='center',
               verticalalignment='center')
    pylab.xlim((1, 1e9))
    pylab.xticks([1, 1e3, 1e6, 1e9])
    pylab.yticks([0, 0.5, 1.0])
    pylab.tight_layout()
    html_writer.embed_matplotlib_figure(fig, width=125, height=125, name="reversibility_cdf")

def main():
    html_fname = '../res/reversibility.html'
    logging.info('Writing HTML output to %s', html_fname)
    html_writer = HtmlWriter(html_fname)
    
    # plot the profile graph
    pylab.rcParams['text.usetex'] = False
    pylab.rcParams['legend.fontsize'] = 10
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 14
    pylab.rcParams['lines.linewidth'] = 2
    pylab.rcParams['lines.markersize'] = 6
    pylab.rcParams['figure.figsize'] = [6.0, 6.0]
    pylab.rcParams['figure.dpi'] = 90
    
    estimators = LoadAllEstimators()
    #analyse_reversibility(estimators['hatzi_gc'], 'HatziGC')
    #analyse_reversibility(estimators['PGC'], 'MiloGC_zoom')
    
    reaction_list = Kegg.getInstance().AllReactions()
    #reaction_list = Feist.FromFiles().reactions
    thermo = estimators['PGC']
    
    thermo.c_mid = DEFAULT_CMID
    thermo.T = DEFAULT_T
    thermo.pH = DEFAULT_PH
    thermo.I = DEFAULT_I
    thermo.pMg = DEFAULT_PMG

    compare_reversibility_to_dG0(reaction_list, thermo=thermo,
                                 html_writer=html_writer)
    
    #compare_reversibility_to_dG0(estimators['alberty'], 'Alberty')
    #metacyc_data('meta', 'zoom', estimators['PGC'], max_pathway_length_for_fig=6)
    #metacyc_data('ecoli', 'zoom', estimators['PGC'], max_pathway_length_for_fig=4)
    
    #meta_regulated_rxns_cumul_plots('meta', '_abs', estimators['PGC'])
    #meta_regulated_rxns_cumul_plots('ecoli', '_abs', estimators['PGC'])
    
    #compare_annotations(reaction_list, thermo=thermo, html_writer=html_writer,
    #                    cmap=GetEmptyConcentrationMap())
    #calc_cons_rxns_corr(estimators['PGC'], 'MiloGC')
    
if __name__ == "__main__":
    main()
