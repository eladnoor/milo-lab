#!/usr/bin/python

import logging
import os
import pylab
import sys

from copy import deepcopy
from matplotlib.font_manager import FontProperties
from optparse import OptionParser
from pygibbs import pathway_modelling
from pygibbs.kegg import Kegg
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from toolbox.util import _mkdir
from pygibbs.pseudoisomer import PseudoisomerMap


class PathwayComparison(object):
    
    def __init__(self, db, html_writer, thermodynamics,
                 kegg=None):
        self.db = db
        self.html_writer = html_writer
        self.thermo = thermodynamics
        self.kegg = kegg or Kegg.getInstance()
        self.pathways = {}
    
    def AddPathway(self, name, pathway_data):
        self.pathways[name] = pathway_data

    def GetConditions(self, pathway_data):
        self.thermo.pH = pathway_data.pH or self.thermo.pH
        self.thermo.I = pathway_data.I or self.thermo.I
        self.thermo.T = pathway_data.T or self.thermo.T
        self.thermo.pMg = pathway_data.pMg or self.thermo.pMg
        self.thermo.c_range = pathway_data.c_range or tuple(self.thermo.c_range)
        self.thermo.c_mid = pathway_data.c_mid or self.thermo.c_mid
        
        self.html_writer.write('Parameters:</br>\n')
        condition_list = ['pH = %g' % self.thermo.pH,
                          'Ionic strength = %g M' % self.thermo.I,
                          'pMg = %g' % self.thermo.pMg,
                          'Temperature = %g K' % self.thermo.T,
                          'Concentration range = %g - %g M' % self.thermo.c_range,
                          'Default concentration = %g M' % self.thermo.c_mid]
        self.html_writer.write_ul(condition_list)

    def GetReactions(self, module_name, pathway_data):
        """Reads the list of reactions from the command file.
        
        Args:
            module_name: the name of the module.
            pathway_data: the pathway.PathwayData object.
        
        Returns:
            A 4 tuple (stoichiometric matrix, rids, fluxes, cids).
        """
        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        cid_mapping = pathway_data.cid_mapping
        field_map = pathway_data.field_map
        
        if pathway_data.kegg_module_id is not None:
            mid = pathway_data.kegg_module_id
            S, rids, fluxes, cids = self.kegg.get_module(mid)
            for i, cid in enumerate(list(cids)):
                if cid in cid_mapping:
                    new_cid, coeff = cid_mapping[cid]
                    cids[i] = new_cid
                    S[:, i] *= coeff
            self.html_writer.write('<h3>Module %s</h3>' % module_name)       
        else:
            S, rids, fluxes, cids = self.kegg.parse_explicit_module(field_map, cid_mapping) 
        
        return (S, rids, fluxes, cids)
     
    def WriteReactionsToHtml(self, S, rids, fluxes, cids, show_cids=True):
        self.thermo.pH = 7
        dG0_r = self.thermo.GetTransfromedReactionEnergies(S, cids)
        
        self.html_writer.write("<li>Reactions:</br><ul>\n")
        
        for r in range(S.shape[0]):
            self.html_writer.write('<li><a href=' + self.kegg.rid2link(rids[r]) + '>R%05d ' % rids[r] + '</a>')
            self.html_writer.write('[x%g, &#x394;G<sub>r</sub><sup>0</sup> = %.1f] : ' % (fluxes[r], dG0_r[r, 0]))
            self.html_writer.write(self.kegg.vector_to_hypertext(S[r, :].flat, cids, show_cids=show_cids))
            self.html_writer.write('</li>\n')
        
        v_total = pylab.dot(pylab.matrix(fluxes), S).flat
        dG0_total = pylab.dot(pylab.matrix(fluxes), dG0_r)[0,0]
        self.html_writer.write('<li><b>Total </b>')
        self.html_writer.write('[&#x394;G<sub>r</sub><sup>0</sup> = %.1f kJ/mol] : \n' % dG0_total)
        self.html_writer.write(self.kegg.vector_to_hypertext(v_total, cids, show_cids=show_cids))
        self.html_writer.write("</li></ul></li>\n")
        
        # Write the kegg-formatted pathway.
        #try:
        #   reactions = map(self.kegg.rid2reaction, rids)
        #   kegg_utils.write_kegg_pathway(self.html_writer, reactions, fluxes)
        #except KeyError, e:
        #   logging.warning('Unknown reaction %s', e)
        
    def PlotPerEnzymeCosts(self, max_dGs, costs):
        costs_array = pylab.array(costs)
        
        fig = pylab.figure()
        pylab.title('Per-enzyme costs', figure=fig)
        pylab.ylabel('Cost per flux unit (log10)', figure=fig)
        pylab.xlabel('Reaction coordinate', figure=fig)
        pylab.semilogy(range(len(costs[0])), costs_array.T, '-*', figure=fig)
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480)
    
    def PlotConcentrations(self, max_dGs, concentrations, cids):
        # Headers, axis labels, ticks for cids
        ncompounds, ntrials = concentrations.shape 
        conc_fig = pylab.figure()
        conc_fig.suptitle('Concentrations')
        pylab.hold(True)
        pylab.xscale('log', figure=conc_fig)
        pylab.ylabel('Compound KEGG ID', figure=conc_fig)
        pylab.xlabel('Concentration [M]', figure=conc_fig)
        
        cids = ["C%05d" % cid for cid in cids]
        pylab.yticks(range(ncompounds, 0, -1), cids,
                     fontproperties=FontProperties(size=8))
        
        # Plot range of concentrations per compound
        for i in xrange(ncompounds):
            concs = concentrations[i,:]
            x_range = (concs.min(), concs.max())
            y_val = ncompounds - i
            y_range = (y_val, y_val + 1e-2)
            
            abs_diff = abs(x_range[0] - x_range[1])
            if abs_diff < 1e-9:
                pylab.scatter(x_range[0], y_range[0], s=1, figure=conc_fig)
            else:
                pylab.plot(x_range, y_range, 'r', figure=conc_fig)
            
            # Plot names
            name = self.kegg.cid2name(cids[i])
            pylab.text(x_range[1] * 1.1, ncompounds - i, name,
                       figure=conc_fig, fontsize=6, rotation=0)
        
        pylab.axis([concentrations.min() / 10, concentrations.max() * 10,
                    0, ncompounds+1],
                   figure=conc_fig)
        
        self.html_writer.embed_matplotlib_figure(conc_fig, width=640, height=480)
    
    def CompareMtdf(self, target_mtdf=None):        
        n_pathways = len(self.pathways)
        for i, (name, pathway_data) in enumerate(self.pathways.iteritems()):
            logging.info('Analyzing pathway %s', name)
            self.html_writer.write('<div margin="20px"><div><b>%s</b></div>' % name)
            self.GetConditions(pathway_data)
            S, rids, fluxes, cids = self.GetReactions(name, pathway_data)
            self.WriteReactionsToHtml(S, rids, fluxes, cids, show_cids=False)
            
            # Bounds on concentrations.         
            bounds = [self.thermo.bounds.get(cid, (None, None))
                      for cid in cids]
            
            # All fluxes are forwards
            fluxes = map(abs, fluxes)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            c_mid = self.thermo.c_mid
            c_range = self.thermo.c_range
            
            path = pathway_modelling.Pathway(S, dG0_f)
            
            if target_mtdf is not None:
                dG_f, conc, score = path.FindMtdf_Regularized(
                    c_range, bounds, c_mid,
                    min_mtdf=target_mtdf,
                    max_mtdf=target_mtdf)
            else:
                dG_f, conc, score = path.FindMTDF_OptimizeConcentrations(
                    c_range, bounds, c_mid)
            if score is None:
                logging.error('No MTDF score for %s', name)
                continue
        
            Nr, Nc = S.shape
            profile_fig = pylab.figure()
            profile_fig.hold(True)
            pylab.title('Thermodynamic Profile',
                        figure=profile_fig)
            pylab.ylabel('Cumulative dG [kJ/mol]', figure=profile_fig)
            pylab.xlabel('Reaction KEGG ID', figure=profile_fig)
            pylab.grid(True, figure=profile_fig)
            
            rids = ['R%05d' % rids[i] for i in xrange(Nr)]
            pylab.xticks(pylab.arange(1, Nr + 1), rids,
                         fontproperties=FontProperties(size=8),
                         rotation=30)
            dG0_r = pylab.zeros((Nr, 1))
            for r in range(Nr):
                reactants = pylab.find(S[r,:])
                dG0_r[r, 0] = pylab.dot(S[r, reactants], dG0_f[reactants])
        
            nan_indices = pylab.find(pylab.isnan(dG0_r))
            finite_indices = pylab.find(pylab.isfinite(dG0_r))
            if (len(nan_indices) > 0):
                dG0_r_finite = pylab.zeros((Nr, 1))
                dG0_r_finite[finite_indices] = dG0_r[finite_indices]
                cum_dG0_r = pylab.cumsum([0] + [dG0_r_finite[r, 0] * fluxes[r] for r in range(Nr)])
            else:
                cum_dG0_r = pylab.cumsum([0] + [dG0_r[r, 0] * fluxes[r] for r in range(Nr)])
            pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG0_r, 'g--', label='Standard [1M]', figure=profile_fig)

            # plot the thermodynamic profile for the different optimization schemes
            dG_r = pylab.dot(S, dG_f)
            self.html_writer.write('<ol>')
            for i, dG in enumerate(dG_r):
                self.html_writer.write('<li>%s: %.2f' % (rids[i], dG))
            self.html_writer.write('</ol>')
            
            cum_dG_r = pylab.cumsum([0] + [dG_r[i, 0] * fluxes[i] for i in range(Nr)])
            pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG_r, figure=profile_fig, label='%s MTDF = %.1f' % (name, score))
            
            pylab.legend(['Standard conditions', 'MTDF'], 'lower left')
            fname = '%s-profile-fig' % name
            
            html_writer.embed_matplotlib_figure(profile_fig, width=640, height=480,
                                                name=fname)

            # Give all compounds with unknown dG0_f the middle concentration value
            conc[nan_indices] = self.thermo.c_mid
            
            unconstrained_cs = []
            unconstrained_cids = []
            for i, bound in enumerate(bounds):
                b_low, b_up = bound
                if b_low is None and b_up is None:
                    unconstrained_cs.append(conc[i, 0])
                    unconstrained_cids.append(cids[i])
            
            n_constrained = len(unconstrained_cs)
            conc_fig = pylab.figure()
            conc_fig.suptitle('Concentrations %s (MTDF = %.1f)' % (name, score))
            pylab.xscale('log', figure=conc_fig)
            pylab.ylabel('Compound KEGG ID', figure=conc_fig)
            pylab.xlabel('Concentration [M]', figure=conc_fig)
            cids_names = ["C%05d" % cid for cid in unconstrained_cids]
            pylab.yticks(range(n_constrained, 0, -1), cids_names,
                         fontproperties=FontProperties(size=8))
            pylab.plot(unconstrained_cs, range(n_constrained, 0, -1),
                       '*b', figure=conc_fig)
    
            x_min = self.thermo.c_range[0] / 10
            x_max = self.thermo.c_range[1] * 50
            y_min = 0
            y_max = n_constrained + 1
            
            for i, concentration in enumerate(unconstrained_cs):
                pylab.text(concentration * 1.1, n_constrained - i,
                           kegg.cid2name(unconstrained_cids[i]),
                           figure=conc_fig, fontsize=6, rotation=0)
                y_val = n_constrained - i
                pylab.plot([x_min, x_max], [y_val, y_val], '-k', linewidth=0.4)
    
            pylab.axvspan(min(unconstrained_cs), max(unconstrained_cs),
                          facecolor='g', alpha=0.3, figure=conc_fig)
            pylab.axis([x_min, x_max, y_min, y_max], figure=conc_fig)
            
            fname = '%s-mtdf-conc-fig' % name
            html_writer.embed_matplotlib_figure(conc_fig, width=640, height=480,
                                                name=fname)

            self.html_writer.write('</div>')

    def CompareMinimumFeasibleConcentrations(self):
        # plot the profile graph
        burdens = {}
        max_dGs = pylab.arange(0.0,-12.25,-0.25)
        
        for name, pathway_data in self.pathways.iteritems():
            self.html_writer.write('<div margin="20px"><div><b>%s</b></div>' % name)
            self.GetConditions(pathway_data)
            S, rids, fluxes, cids = self.GetReactions(name, pathway_data)
            self.WriteReactionsToHtml(S, rids, fluxes, cids, show_cids=False)
            
            # Bounds on concentrations.            
            bounds = [self.thermo.bounds.get(cid, (None, None))
                      for cid in cids]
            
            # All fluxes are forwards
            fluxes = map(abs, fluxes)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            c_mid = self.thermo.c_mid
            c_range = self.thermo.c_range
            
            path = pathway_modelling.Pathway(S, dG0_f)
            dgs, concentrations, total_conc = path.FindKineticOptimum()
            
            self.html_writer.write('<div margin="20px"><div>Total concentration <b>%f M</b></div>' % total_conc)
                            
            concs_array = pylab.array(concentrations)
            self.PlotConcentrations(max_dGs, concs_array, cids)
            
            self.html_writer.write('</div>')
        
        
    def CompareMinimalEnzymeBurden(self):
        # plot the profile graph
        burdens = {}
        max_dGs = pylab.arange(0.0,-12.25,-0.25)
        
        for name, pathway_data in self.pathways.iteritems():
            self.html_writer.write('<div margin="20px"><div><b>%s</b></div>' % name)
            self.GetConditions(pathway_data)
            S, rids, fluxes, cids = self.GetReactions(name, pathway_data)
            self.WriteReactionsToHtml(S, rids, fluxes, cids, show_cids=False)
            
            # Bounds on concentrations.            
            bounds = [self.thermo.bounds.get(cid, (None, None))
                      for cid in cids]
            
            # All fluxes are forwards
            fluxes = map(abs, fluxes)
            dG0_f = self.thermo.GetTransformedFormationEnergies(cids)
            c_mid = self.thermo.c_mid
            c_range = self.thermo.c_range
            
            path = pathway_modelling.Pathway(S, dG0_f)
            costs = []
            per_enzyme_costs = []
            concentrations = []
            for max_dg in max_dGs:
                total_cost, per_enzyme, concs = path.FindPcrEnzymeCost(
                    c_mid=c_mid, ratio=3.0, bounds=bounds, max_reaction_dg=max_dg,
                    fluxes=fluxes)
                costs.append(total_cost)
                per_enzyme_costs.append(per_enzyme)
                concentrations.append(concs)
                
            burdens[name] = costs
            concs_array = pylab.hstack(concentrations)
            self.PlotPerEnzymeCosts(max_dGs, per_enzyme_costs)
            self.PlotConcentrations(max_dGs, concs_array, cids)
            
            
            self.html_writer.write('</div>')
        
        self.html_writer.write('<div margin="20px"><div><b>Minimal Enzyme Costs</b></div>')
        self.html_writer.write('<table border="1px"><tr><th>Max dG</th>')
        for dG in max_dGs:
            self.html_writer.write('<th>%.3f</th>' % dG)
        self.html_writer.write('</tr>')
        
        for name, costs in burdens.iteritems():
            self.html_writer.write('<tr>')
            self.html_writer.write('<td>%s</td>' % name)
            for cost in costs:
                self.html_writer.write('<td>%.4f</td>' % cost)
            self.html_writer.write('</tr>')
        self.html_writer.write('</table>')
        self.html_writer.write('</div>')
        
        fig = pylab.figure()
        pylab.hold(True)
        keys = burdens.keys()
        pathway_datas = map(self.pathways.get, keys)
        titles = [p.name for p in pathway_datas]
        minus_max_dGs = -pylab.array(max_dGs)
        costs = pylab.array(map(burdens.get, keys))
        pylab.semilogy(minus_max_dGs, costs.T, '-', figure=fig)
        pylab.axvline(1.0, linestyle='--', color='k', linewidth=2, figure=fig)
        pylab.axvline(2.7, linestyle='--', color='k', linewidth=2, figure=fig)
        pylab.axvline(5.7, linestyle='--', color='k', linewidth=2, figure=fig)
        pylab.axvline(11.4, linestyle='--', color='k', linewidth=2, figure=fig)
        
        pylab.legend(titles, loc='upper left')
        pylab.ylabel("Burden per input flux unit (log10)", figure=fig)
        pylab.xlabel("-Max dG [kJ/mol]", figure=fig)
        
        name = 'enzyme-burden'
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480,
                                                 name=name)



def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=['observed_only',
                                   'hatzi_only',
                                   'milo_only',
                                   'milo_merged'],
                          default="milo_merged",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermodynamics_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="The name of the thermodynamics file to load.")
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          default="../data/thermodynamics/pathways.txt",
                          help="The file to read for pathways to analyze.")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          default='../res/thermo_comparison/report.html',
                          help="Where to write output to.")
    return opt_parser


if __name__ == "__main__":    
    options, _ = MakeOpts().parse_args(sys.argv)
    input_filename = os.path.abspath(options.input_filename)
    output_filename = os.path.abspath(options.output_filename)
    if not os.path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename
    print 'Will write output to %s' % output_filename
    
    thermo_source = options.thermodynamics_source
    assert thermo_source in ('observed_only', 'hatzi_only', 'milo_only', 'milo_merged')
    print 'Will load thermodynamic data from source "%s"' % thermo_source
    
    db_loc = options.db_filename
    print 'Reading from DB %s' % db_loc
    db = SqliteDatabase(db_loc)
    
    public_db_loc = options.kegg_db_filename
    print 'Reading from public DB %s' % public_db_loc
    db_public = SqliteDatabase(public_db_loc)
    
    # dG0 =  -E'*nE*F - R*T*ln(10)*nH*pH
    # Where: 
    #    F  = 96.48 kC/mol
    #    nE - change in e-
    #    nH - change in H+
    #    pH - the conditions in which the E' was measured
    #
    # Ferredoxin  ox/red: E' = -0.380V (nE = 1, nH = 0) -> dG0 = 38.0 kJ/mol [1]
    # Ubiqinone   ox/red: E' =  0.113V (nE = 2, nH = 2) -> dG0 = -103.2 kJ/mol [1]
    # Menaquinone ox/red: E' = -0.074V (nE = 2, nH = 2) -> dG0 = -65.8 kJ/mol [1]
    #
    # [1] - Thauer 1977
    
    observed_thermo_fname = options.thermodynamics_filename
    print 'Loading observed thermodynamic data from %s' % observed_thermo_fname
    observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
        observed_thermo_fname)
    
    if thermo_source == 'hatzi_only':    
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(db, 'hatzi_thermodynamics')
        thermo.AddPseudoisomer( 139, nH=0,  z=1, nMg=0, dG0=0)      # Ferrodoxin(ox)
        thermo.AddPseudoisomer( 138, nH=0,  z=0, nMg=0, dG0=38.0)   # Ferrodoxin(red)
        thermo.AddPseudoisomer( 399, nH=90, z=0, nMg=0, dG0=0)      # Ubiquinone-10(ox)
        thermo.AddPseudoisomer( 390, nH=92, z=0, nMg=0, dG0=-103.2) # Ubiquinone-10(red)
        thermo.AddPseudoisomer( 828, nH=16, z=0, nMg=0, dG0=0)      # Menaquinone(ox)
        thermo.AddPseudoisomer(5819, nH=18, z=0, nMg=0, dG0=-65.8)  # Menaquinone(red)
        thermo.SetPseudoisomerMap(101, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=0.0)) # THF
        thermo.SetPseudoisomerMap(234, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=-137.5)) # 10-Formyl-THF
        thermo.SetPseudoisomerMap(445, PseudoisomerMap(nH=22, z=0, nMg=0, dG0=65.1)) # 5,10-Methenyl-THF
        thermo.SetPseudoisomerMap(143, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=77.9)) # 5,10-Methylene-THF
        thermo.SetPseudoisomerMap(440, PseudoisomerMap(nH=25, z=0, nMg=0, dG0=32.1)) # 5-Methyl-THF
    elif thermo_source == 'milo_only':
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(
            db, 'pgc_pseudoisomers')
    elif thermo_source == 'milo_merged':
        thermo = PsuedoisomerTableThermodynamics.FromDatabase(
            db, 'pgc_pseudoisomers')
        thermo.override_data(observed_thermo)
    elif thermo_source == 'observed_only':
        thermo = observed_thermo
    else:
        logging.fatal('Unknown thermodynamic data source.')
    
    kegg = Kegg.getInstance()
    thermo.bounds = deepcopy(kegg.cid2bounds)
    thermo.c_range = (1e-7,1e-1)
    
    dirname = os.path.dirname(output_filename)
    if not os.path.exists(dirname):
        print 'Making output directory %s' % dirname
        _mkdir(dirname)
    
    print 'Executing thermodynamic pathway analysis'
    html_writer = HtmlWriter(output_filename)
    comparator = PathwayComparison(db, html_writer, thermodynamics=thermo)
    
    entry2fields_map = ParsedKeggFile.FromKeggFile(input_filename)
    for name, field_map in entry2fields_map.iteritems():
        pathway_data = PathwayData.FromFieldMap(field_map)
        if pathway_data.skip:
            logging.info('Skipping pathway %s', pathway_data.name)
            continue
        
        comparator.AddPathway(name, pathway_data)
        
    comparator.CompareMtdf()
    
