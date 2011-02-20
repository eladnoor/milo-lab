from scipy import stats
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs import thermodynamics
from pygibbs.thermodynamic_constants import R, default_I, default_pH
from pygibbs.thermodynamic_constants import default_pMg, default_T
from pygibbs.groups import GroupContribution
import pylab
from pygibbs.kegg import Kegg,KeggNonCompoundException
from toolbox.plotting import cdf
from SOAPpy import WSDL 
from pygibbs.metacyc import MetaCyc,MetaCycNonCompoundException
import logging
import matplotlib

def try_kegg_api():
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/dG0_test.html')
    kegg = Kegg(db)
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
            sparse = G.kegg().rid2sparse_reaction(rid)
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
    return cmap


def ConcentrationFactor(sparse_reaction,
                        concentration_map,
                        c_mid):
    factor = 0
    for cid, stoic in sparse_reaction.iteritems():
        concentration = concentration_map.get(cid, None) or c_mid
        factor += pylab.log(concentration) * stoic
    return factor

def CalculateReversability(sparse, G, c_mid=1e-3, pH=default_pH, 
                           pMg=default_pMg, I=default_I, T=default_T,
                           concentration_map=None):
    cmap = concentration_map or {}
    dG0 = G.estimate_dG_reaction(sparse, pH, pMg, I, T)
    
    # remove H2O and H+ from the list of reactants since their
    # concentration is fixed
    sparse.pop(WATER, None)
    sparse.pop(HPLUS, None)

    cfactor = ConcentrationFactor(sparse, cmap, c_mid)
    sum_abs_s = sum([abs(x) for k, x in sparse.iteritems()
                     if k not in cmap])
    if (sum_abs_s == 0):
        return None
    else:
        return 2 / pylab.log(10) * ((-dG0/(R*T) + cfactor) / sum_abs_s)

def calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, cmap):
    histogram = {}
    histogram['total'] = []
    hits = 0
    misses = 0
    debug_file = open('../res/kegg_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Module\tPosition\tReaction Name\tDefinition\tEC list\tEquation\tRev IND\n")
    
    for mid, rid_flux_list in kegg.mid2rid_map.iteritems():
        if not rid_flux_list or len(rid_flux_list) < 2:
            continue
        for i, (rid, flux) in enumerate(rid_flux_list):
            try:
                r = flux * CalculateReversability(G.kegg().rid2sparse_reaction(rid), G, c_mid, pH, pMg, I, T,
                                                  concentration_map=cmap)
                rxn = kegg.rid2reaction_map[rid]
                debug_file.write("%s\t%d\t%s\t%s\t%s\t%s\t%f\n" % (kegg.mid2name_map[mid], i+1, rxn.name, rxn.definition, rxn.ec_list, rxn.equation, r))
                histogram.setdefault(i, []).append(r)
                if i > 1:
                    histogram['total'].append(r)
                hits += 1
            except thermodynamics.MissingCompoundFormationEnergy:
                misses += 1
                continue
    debug_file.close()
    
    logging.info("Reactions with known dG0: %d" % hits)
    logging.info("Reactions with unknown dG0: %d" % misses)
    return histogram

def plot_histogram(histogram, html_writer, title='', max_pathway_length=8):
    fig = pylab.figure()
    pylab.hold(True)

    colors = {0:'r', 1:'orange', 2:'green', 3:'cyan', 4:'blue', 5:'violet', 
              6:'pink', 'total':'k--'}
    for key, value in histogram.iteritems():
        if len(value) > 15:
            cdf(value, label='%s (median=%.1f, N=%d)' % \
                (key, stats.cmedian(value), len(value)),
                style=colors.get(key, 'grey'))
    pylab.xlim(-20, 20)
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

def calculate_metacyc_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, metacyc, cmap):
    histogram = {}
    histogram['total'] = []
    hits = 0
    misses = 0
    comp_translation_misses = 0
    flux_misses = 0
    rxn_map_misses = 0
    meta_compound_misses = 0
    
    debug_file = open('../res/metacyc_' + ('constrained_' if len(cmap) > 0 else 'non_constrained_') + 'rev.txt', 'w')
    debug_file.write("Pathway\tPosition\tReaction Name\tEC list\tEquation\tRev IND\n")
    
    for pathway in metacyc.uid2pathway_map.itervalues():
        rxns_dict = pathway.GetRxnsOrder()
        rxn_dirs = pathway.GetRxnsDirs()
        
        if (len(rxns_dict) >= 2):
            for rxn,pos in rxns_dict.iteritems():
                try:
                    sparse = metacyc.rxn_uid2sparse_reaction(rxn)
                    if (not sparse):
                        rxn_map_misses += 1
                        continue
                    if (rxn in rxn_dirs):
                        flux = rxn_dirs[rxn]
                    else:
                        flux_misses += 1
                        continue
                    
                    r = CalculateReversability(metacyc.sparse2kegg_cids(sparse, kegg), G, c_mid, pH, pMg, I, T,
                                                      concentration_map=cmap)
                    if (r != None):
                        r *= flux
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
                
    debug_file.close()
    logging.info("Reactions with known dG0: %d" % hits)
    logging.info("Reactions with unknown dG0: %d" % misses)
    logging.info("Reactions with unknown compounds: %d" % comp_translation_misses)
    logging.info("Reactions with missing flux direction: %d" % flux_misses)
    logging.info("Reactions not in MetaCyc reactions file: %d" % rxn_map_misses)
    logging.info("Reactions with compounds not in compounds file: %d" % meta_compound_misses)

    return histogram

def main():
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/reversibility.html')
    kegg = Kegg(db)
    G = GroupContribution(db, html_writer=html_writer, kegg=kegg)
    G.init()
    c_mid = 1e-3
    pH, pMg, I, T = (7.0, 3.0, 0.1, 298.15)
    
    histogram = calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg,
                                                  cmap=GetConcentrationMap(kegg))
    fig1 = plot_histogram(histogram, html_writer, title='With constraints on co-factors')
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    
    histogram = calculate_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg,
                                                  cmap={})
    fig2 = plot_histogram(histogram, html_writer, title='No constraints on co-factors')
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)

def metacyc_data(org):
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/' + org + '_reversibility.html')
    kegg = Kegg(db)
    metacyc_inst = MetaCyc(org, db)
    G = GroupContribution(db, html_writer=html_writer, kegg=kegg)
    G.init()
    c_mid = 1e-3
    pH, pMg, I, T = (7.0, 3.0, 0.1, 298.15)
    
    histogram = calculate_metacyc_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, metacyc_inst,
                                                  cmap=GetConcentrationMap(kegg))
    fig1 = plot_histogram(histogram, html_writer, title=('%s pathways: With constraints on co-factors' % org))
    html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility1.pdf', figure=fig1, format='pdf')
    
    histogram = calculate_metacyc_reversibility_histogram(G, c_mid, pH, pMg, I, T, kegg, metacyc_inst,
                                                  cmap={})
    fig2 = plot_histogram(histogram, html_writer, title=('%s pathways: No constraints on co-factors' % org ))
    html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
    pylab.savefig('../res/' + org + '_reversibility2.pdf', figure=fig1, format='pdf')
    
if __name__ == "__main__":
    #try_kegg_api()
    main()
    metacyc_data('meta')
    #metacyc_data('ecoli')