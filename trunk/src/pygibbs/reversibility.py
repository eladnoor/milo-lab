#!/usr/bin/python

from scipy import stats
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs import thermodynamics
from pygibbs.thermodynamic_constants import R, default_I, default_pH
from pygibbs.thermodynamic_constants import default_pMg, default_T
from pygibbs.groups import GroupContribution
import pylab
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggNonCompoundException,\
    KeggReactionNotBalancedException
from toolbox.plotting import cdf
from SOAPpy import WSDL 
from pygibbs.metacyc import MetaCyc,MetaCycNonCompoundException
import logging
import matplotlib
import pygibbs.kegg_utils
import re
import math
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
import numpy

def try_kegg_api():
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/dG0_test.html')
    kegg = Kegg.getInstance()
    G = GroupContribution(db, html_writer=html_writer, kegg=kegg)
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
    cdf(rid2reversibility.values(), 'all reactions', 'r', show_median=True)
    pylab.show()
    

WATER = 1
HPLUS = 80


def GetConcentrationMap(kegg_handle):
    cmap = {}    
    for cid in kegg_handle.get_all_cids():
        lower, upper = kegg_handle.get_bounds(cid)
        if lower and upper:
            # In the file we got this data from lower = upper 
            cmap[cid] = lower
    cmap[WATER] = 1
    cmap[HPLUS] = 1
    return cmap


def ConcentrationFactor(sparse_reaction,
                        concentration_map,
                        c_mid):
    factor = 0
    for cid, stoic in sparse_reaction.iteritems():
        concentration = concentration_map.get(cid, None) or c_mid
        factor += pylab.log(concentration) * stoic
    return factor

def CalculateReversability(sparse, thermo, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T,
                           concentration_map=None):
    cmap = concentration_map or {}
    dG0 = thermo.reaction_to_dG0(sparse, pH, pMg, I, T)
    
    cfactor = ConcentrationFactor(sparse, cmap, c_mid)
    sum_abs_s = sum([abs(x) for k, x in sparse.iteritems()
                     if k not in cmap])
    if (sum_abs_s == 0):
        return None
    else:
        return 2 / pylab.log(10) * ((-dG0/(R*T) - cfactor) / sum_abs_s)

def calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, cmap, id):
    histogram = {}
    histogram['total'] = []
    
    rel_histogram = {}
    rel_histogram['total'] = []
    
    n_first_max = 0
    n_valid = 0
    
    hits = 0
    misses = 0
    non_balanced = 0
    
    debug_file = open('../res/kegg_' + id + '_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Module\tPosition\tReaction Name\tDefinition\tEC list\tEquation\tRev IND\n")
    
    for mid, rid_flux_list in kegg.mid2rid_map.iteritems():
        if not rid_flux_list or len(rid_flux_list) < 2:
            continue
        
        pw_r_map = {}
        pw_r_sum = 0
        pw_n_r = 0
        
        for i, (rid, flux) in enumerate(rid_flux_list):
            try:
                sparse = kegg.rid2sparse_reaction(rid)
                sparse = kegg.BalanceReaction(sparse, balance_water=True)
                r = flux * CalculateReversability(sparse, G, c_mid, pH, pMg, I, T,
                                                  concentration_map=cmap)
                rxn = kegg.rid2reaction_map[rid]
                dbg = "%s__DEL__%d__DEL__%s__DEL__%s__DEL__%s__DEL__%s__DEL__%f\n" % (kegg.mid2name_map[mid], i+1, rxn.name, rxn.definition, rxn.ec_list, rxn.equation, r)
                dbg = re.sub ('\t', ' ', dbg)
                debug_file.write(re.sub('__DEL__', '\t', dbg))
                histogram.setdefault(i, []).append(r)
                if i > 1:
                    histogram['total'].append(r)
                hits += 1
                
                pw_n_r += 1
                pw_r_sum += abs(r)
                if (i in pw_r_map):
                    pw_r_map[i] = pw_r_map[i] + ';' + str(r)
                else:
                    pw_r_map[i] = str(r)
                            
            except thermodynamics.MissingCompoundFormationEnergy:
                misses += 1
                continue
            except KeggReactionNotBalancedException:
                    non_balanced += 1
                    #print 'Reaction cannot be balanced, uid: %s' % rxn
                    continue
            
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
                        if pos > 0:
                            rel_histogram['total'].append(r)
                if (first_r == curr_max):
                        n_first_max += 1
                                    
    debug_file.close()
    
    logging.info("Reactions with known dG0: %d" % hits)
    logging.info("Reactions with unknown dG0: %d" % misses)
    logging.info("Non balanced reactions: %d" % non_balanced)
    
    return histogram,rel_histogram,(float(n_first_max)/n_valid)

def plot_histogram(histogram, html_writer, title='', max_pathway_length=8, xlim=20):
    fig = pylab.figure()
    pylab.hold(True)

    colors = {0:'r', 1:'orange', 2:'green', 3:'cyan', 4:'blue', 5:'violet', 
              6:'pink', 'total':'k--'}
    for key, value in histogram.iteritems():
        if len(value) > 15:
            cdf(value, label='%s (median=%.1f, N=%d)' % \
                (key, stats.cmedian(value), len(value)),
                style=colors.get(key, 'grey'))
    pylab.xlim(-1 * xlim, xlim)
    pylab.xlabel('irreversability')
    pylab.ylabel('cumulative distribution')
    legendfont = matplotlib.font_manager.FontProperties(size=7)
    pylab.legend(loc='lower right', prop=legendfont)
    pylab.title(title)
    pylab.hold(False)
    
    print stats.ranksums(histogram[0], histogram['total'])
    #for k1, h1 in histogram.iteritems():
    #    for k2, h2 in histogram.iteritems():
    #        print k1, k2, stats.ranksums(h1, h2)
    
    return fig


def calculate_metacyc_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, kegg, metacyc, cmap, org):
    histogram = {}
    histogram['total'] = []
    
    rel_histogram = {}
    rel_histogram['total'] = []
    
    n_first_max = 0
    n_valid = 0
    
    hits = 0
    misses = 0
    comp_translation_misses = 0
    flux_misses = 0
    rxn_map_misses = 0
    meta_compound_misses = 0
    meta_non_balanced = 0
    
    debug_file = open('../res/metacyc_' + org + ('_constrained_' if len(cmap) > 0 else '_non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Pathway\tPosition\tReaction Name\tEC list\tEquation\tRev IND\n")
    
    for pathway in metacyc.uid2pathway_map.itervalues():
        rxns_dict = pathway.GetRxnsOrder()
        rxn_dirs = pathway.GetRxnsDirs()
        
        curr_max = -1000
        
        pw_r_map = {}
        pw_r_sum = 0
        pw_n_r = 0
        
        if (len(rxns_dict) >= 2):
            for rxn,pos in rxns_dict.iteritems():
                try:
                    sparse = metacyc.rxn_uid2sparse_reaction(rxn)
                    
                    if (not sparse):
                        #print 'Reaction not found, uid: %s' % rxn
                        rxn_map_misses += 1
                        continue
                    
                    sparse = kegg.BalanceReaction(metacyc.sparse2kegg_cids(sparse, kegg), balance_water=True)
                    
                    if (rxn in rxn_dirs):
                        flux = rxn_dirs[rxn]
                    else:
                        flux_misses += 1
                        continue
                    
                    r = CalculateReversability(sparse, thermo, c_mid, pH, pMg, I, T,
                                                      concentration_map=cmap)
                    if (r != None):
                        r *= flux
                                                
                        pw_n_r += 1
                        pw_r_sum += abs(r)
                        if ((pos - 1) in pw_r_map):
                            pw_r_map[pos - 1] = pw_r_map[pos - 1] + ';' + str(r)
                        else:
                            pw_r_map[pos - 1] = str(r)
                            
                        debug_file.write("%s\t%d\t%s\t%s\t%s\t%f\n" % (str(pathway.name), pos, metacyc.uid2reaction_map[rxn].name, metacyc.uid2reaction_map[rxn].ec_number, metacyc.uid2reaction_map[rxn].equation, r))
                        histogram.setdefault(pos - 1, []).append(r)
                        if pos >= 1:
                            histogram['total'].append(r)
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
                            if pos > 0:
                                rel_histogram['total'].append(r)
                    if (first_r == curr_max):
                        n_first_max += 1
                
    debug_file.close()
    logging.info("Reactions with known dG0: %d" % hits)
    logging.info("Reactions with unknown dG0: %d" % misses)
    logging.info("Reactions with unknown compounds (translation to Kegg failed): %d" % comp_translation_misses)
    logging.info("Reactions with missing flux direction: %d" % flux_misses)
    logging.info("Reactions with compounds not found in MetaCyc compounds file: %d" % rxn_map_misses)
    logging.info("Reactions with compounds not in MetaCyc compounds file: %d" % meta_compound_misses)
    logging.info("Reactions which couldn't be balanced: %d" % meta_non_balanced)

    return histogram,rel_histogram,(float(n_first_max)/n_valid)

def main(thermo, name):
    html_writer = HtmlWriter('../res/' + name + '_reversibility.html')
    kegg = Kegg.getInstance()
    c_mid = 1e-3
    pH, pMg, I, T = (7.0, 3.0, 0.1, 298.15)
    
    (histogram,rel_histogram,perc_first_max) = calculate_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, kegg,
                                                  cmap=GetConcentrationMap(kegg), id=name)
    
    html_writer.write('<h1>' + name + ': Constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig1 = plot_histogram(histogram, html_writer, title='With constraints on co-factors', xlim=20)
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility1.png', figure=fig1, format='png')
    
    fig1_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module with constraints on co-factors', xlim=5)
    html_writer.embed_matplotlib_figure(fig1_rel, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility1_rel.png', figure=fig1_rel, format='png')
    
    (histogram,rel_histogram, perc_first_max) = calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg,
                                                  cmap={}, id=name)

    html_writer.write('<h1>' + name + ': Non constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig2 = plot_histogram(histogram, html_writer, title='No constraints on co-factors', xlim=20)
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2.png', figure=fig2, format='png')
    
    fig2_rel = plot_histogram(rel_histogram, html_writer, title='Normed per module, no constraints on co-factors', xlim=5)
    html_writer.embed_matplotlib_figure(fig2_rel, width=640, height=480)
    pylab.savefig('../res/' + name + '_kegg_reversibility2_rel.png', figure=fig2_rel, format='png')
    
def metacyc_data(org, thermo):
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/' + org + '_reversibility.html')
    kegg = Kegg.getInstance()
    metacyc_inst = MetaCyc(org, db)
    c_mid = 1e-3
    pH, pMg, I, T = (7.0, 3.0, 0.1, 298.15)
    
    (histogram,rel_histogram,perc_first_max) = calculate_metacyc_reversibility_histogram(thermo, c_mid, pH, pMg, I, T, kegg, metacyc_inst,
                                                  cmap=GetConcentrationMap(kegg), org=org)
    
    html_writer.write('<h1>Constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig1 = plot_histogram(histogram, html_writer, title=('%s pathways: With constraints on co-factors' % org), xlim=20)
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility1.png', figure=fig1, format='png')

    fig1_rel = plot_histogram(rel_histogram, html_writer, title=('%s pathways: Normed per pathway with constraints on co-factors' % org), xlim=5)
    html_writer.embed_matplotlib_figure(fig1_rel, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility1_rel.png', figure=fig1_rel, format='png')
    
    (histogram,rel_histogram,perc_first_max) = calculate_metacyc_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, metacyc_inst,
                                                  cmap={}, org=org)
    
    html_writer.write('<h1>Non constrained co-factors</h1>Percentage of modules where first reaction is the maximal: %f<br>' % perc_first_max)
    fig2 = plot_histogram(histogram, html_writer, title=('%s pathways: No constraints on co-factors' % org ), xlim=20)
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility2.png', figure=fig1, format='png')

    fig2_rel = plot_histogram(rel_histogram, html_writer, title=('%s pathways: Normed per pathway no constraints on co-factors' % org), xlim=5)
    html_writer.embed_matplotlib_figure(fig2_rel, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility2_rel.png', figure=fig2_rel, format='png')

def calc_palsson_rxns_irrev(org, rxns_file, comp2cid_file, thermo, c_mid, pH, pMg, I, T, cmap):

    html_writer = HtmlWriter('../res/' + org + '_model_rev.html')
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
                n_hits += 1
                histogram.setdefault(reversible, []).append(r)
            
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
        
    
    logging.info("Reactions with known dG0: %d" % n_hits)
    logging.info("Reactions with unknown dG0: %d" % misses)
    logging.info("Reactions with unknown compounds (translation to Kegg failed): %d" % n_miss_comp)
    logging.info("Non balanced reactions: %d" % n_not_balanced)
    logging.info("Reactions of solely currency copmounds: %d" % n_only_currency)
    
    # plot the bar 
    fig = pylab.figure()
    pylab.hold(True)

    colors = {'Reversible':'blue', 'Irreversible':'red'}
    
    i = 1
    for key, value in histogram.iteritems():    
        pylab.subplot (2,1,i)
        i +=1
        label = '%s (%d reactions)' % (key, len(value))      
        pylab.hist(value, 101, range=(-50,50), normed=True, label=label, color=colors[key])
        legendfont = matplotlib.font_manager.FontProperties(size=9)
        pylab.legend(prop=legendfont)
        pylab.xlabel('irreversability')
        pylab.ylabel('Fraction of reactions')
        #pylab.title(org + ' model')
    
    pylab.hold(False)
    
    html_writer.write('<h1>Irreversibility index histogram on Palsson\'s model %s</h1>' % org)
    html_writer.embed_matplotlib_figure(fig, width=640, height=480)
    pylab.savefig('../res/' + org + '_model_irrev.png', figure=fig, format='png')
    
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
            sparse[int(cid[1:])] = coef * s    
        else:
            raise KeggNonCompoundException
        
class IrrevParseException(Exception):
    pass

if __name__ == "__main__":
    #logging.getLogger('').setLevel(logging.DEBUG)
    #try_kegg_api()

    db_public = SqliteDatabase('../res/gibbs.sqlite')
    G = GroupContribution(db_public)
    G.init()
    kegg = Kegg.getInstance()
    G_hatzi = Hatzi()
    #G = PsuedoisomerTableThermodynamics.FromDatabase(db_public, 'alberty_pseudoisomers')
    
    main(G, 'MiloGC')
    #main(G_hatzi, 'HatziGC')
    metacyc_data('meta', G)
    metacyc_data('ecoli', G)
    
    
    #calc_palsson_rxns_irrev('ecoli', '../data/reactionList_iAF1260.txt', '../data/ecoli2kegg_cid.txt', thermo=G, c_mid=1e-3, pH=7.0, pMg=3.0, I=0.1, T=298.15, cmap=GetConcentrationMap(kegg))