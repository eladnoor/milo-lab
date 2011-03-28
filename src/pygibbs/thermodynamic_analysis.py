import re
from thermodynamic_constants import default_T, default_pH, default_I,\
    default_c0, R
import pylab
import copy
from pygibbs.feasibility import pC_to_range, find_mtdf, find_pCr,\
    LinProgNoSolutionException, thermodynamic_pathway_analysis
from pygibbs.thermodynamic_constants import transform
import matplotlib
import logging
from toolbox.html_writer import HtmlWriter
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.kegg import Kegg
from pygibbs import kegg_utils

class ThermodynamicAnalysis(object):
    def __init__(self, db, html_writer, thermodynamics):
        self.thermo = thermodynamics
        self.html_writer = html_writer
        self.kegg = Kegg.getInstance()
        self.db = db

    def analyze_pathway(self, filename):
        self.html_writer.write("<h1>Pathway analysis using Group Contribution Method</h1>\n")
        entry2fields_map = ParsedKeggFile.FromKeggFile(filename)
        
        for key in sorted(entry2fields_map.keys()):
            field_map = entry2fields_map[key]
            if field_map.GetBoolField('SKIP'):
                logging.info("skipping pathway: " + key)
                continue
            try:
                self.html_writer.write("<b>%s - %s</b>" % (field_map.GetStringField('NAME'), field_map.GetStringField('TYPE')))
                self.html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (key))
                self.html_writer.write('<div id="%s" style="display:none">' % key)
                self.html_writer.write('<h2>%s - %s</h2>\n' % (field_map.GetStringField('NAME'), field_map.GetStringField('TYPE')))
            except KeyError:
                raise Exception("Both the 'NAME' and 'TYPE' fields must be defined for each pathway")

            logging.info("analyzing pathway: " + key)

            function_dict = {'PROFILE':self.analyze_profile,
                             'SLACK':self.analyze_slack,
                             'MARGIN':self.analyze_margin,
                             'CONTOUR':self.analyze_contour,
                             'PROTONATION':self.analyze_protonation}

            analysis_type = field_map.GetStringField('TYPE')
            if analysis_type in function_dict:
                function_dict[analysis_type](key, field_map)     
            else:
                raise Exception("Unknown analysis type: " + analysis_type)
            
            self.html_writer.write('</div><br>\n')
            
        self.html_writer.write('<h4>Measured concentration table:</h4>\n')
        self.html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % ('__concentrations__'))
        self.html_writer.write('<div id="%s" style="display:none">' % '__concentrations__')
        self.html_writer.write('<p><h2>Abundance</h2>\n')
        self.db.Query2HTML(self.html_writer,
                           "SELECT cid, media, 1000*concentration from compound_abundance ORDER BY cid, media",
                           column_names=["cid", "media", "concentration [mM]"])
        self.html_writer.write('</p>\n')
        self.html_writer.write('</div><br>\n')

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if (len(tokens) == 0):
            return default_value
        if (len(tokens) > 1):
            raise Exception("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    def get_reactions(self, module_name, field_map):
        """
            read the list of reactions from the command file
        """
        
        if "MODULE" in field_map:
            mid_str = field_map["MODULE"]
            if (mid_str[0] == 'M'):
                mid = int(mid_str[1:])
            else:
                mid = int(mid_str)
            (S, rids, fluxes, cids) = self.kegg.get_module(mid)
            self.html_writer.write('<h3>Module <a href=http://www.genome.jp/dbget-bin/www_bget?M%05d>M%05d</a></h3>\n' % (mid, mid))       
        else:
            (S, rids, fluxes, cids) = self.kegg.parse_explicit_module(field_map) 
        
        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        if "MAP_CID" in field_map:
            for line in field_map["MAP_CID"].split('\t'):
                (cid_before, cid_after) = [int(cid[1:]) for cid in line.split(None, 1)]
                if (cid_before in cids):
                    cids[cids.index(cid_before)] = cid_after
        
        return (S, rids, fluxes, cids)

    def write_metabolic_graph(self, name, S, rids, cids):
        """
            draw a graph representation of the pathway
        """        
        Gdot = self.kegg.draw_pathway(S, rids, cids)
        self.html_writer.embed_dot(Gdot, name, width=400, height=400)

    def analyze_profile(self, key, field_map):
        self.html_writer.write('<p>\n')
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # read the list of conditions from the command file
        conditions = []
        for condition in field_map["CONDITIONS"].split('\t'):
            (media, pH, I, T, c0) = (None, default_pH, default_I, default_T, default_c0)
            media = re.findall("media=([a-zA-Z_]+)", condition)[0]
            if (media == 'None'):
                media = None
            pH = ThermodynamicAnalysis.get_float_parameter(condition, "pH", default_pH)
            I = ThermodynamicAnalysis.get_float_parameter(condition, "I", default_I)
            T = ThermodynamicAnalysis.get_float_parameter(condition, "T", default_T)
            c0 = ThermodynamicAnalysis.get_float_parameter(condition, "c0", default_c0)
            conditions.append((media, pH, I, T, c0))
            self.html_writer.write('<li>Conditions: media = %s, pH = %g, I = %g M, T = %g K, c0 = %g</li>\n' % (media, pH, I, T, c0))
        self.html_writer.write('</ol></li>\n')
        
        # read the list of methods for calculating the dG
        methods = []
        if field_map.GetBoolField('MILO', default_value=True):
            methods.append('MILO')
        if field_map.GetBoolField('HATZI', default_value=False):
            methods.append('HATZI')
        
        # prepare the legend for the profile graph
        legend = []
        dG_profiles = {}
        params_list = []
        for (media, pH, I, T, c0) in conditions:
            for method in methods:
                plot_key = method + ' dG (media=%s,pH=%g,I=%g,T=%g,c0=%g)' % (str(media), pH, I, T, c0)
                legend.append(plot_key)
                dG_profiles[plot_key] = []
                params_list.append((method, media, pH, I, T, c0, plot_key))

        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>')
        self.write_metabolic_graph(key, S, rids, cids)
        
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = {}
        dG0_r = {}
        dG_f = {}
        dG_r = {}
        for (method, media, pH, I, T, c0, plot_key) in params_list:
            dG0_f[plot_key] = pylab.zeros((Nc, 1))
            dG_f[plot_key] = pylab.zeros((Nc, 1))
            for c in range(Nc):
                if method == "MILO":
                    dG0_f[plot_key][c] = self.thermo.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                elif method == "HATZI":
                    dG0_f[plot_key][c] = self.thermo.hatzi.cid2dG0_tag(cids[c], pH=pH, I=I, T=T)
                else:
                    raise Exception("Unknown dG evaluation method: " + method)
                # add the effect of the concentration on the dG_f (from dG0_f to dG_f)
                dG_f[plot_key][c] = dG0_f[plot_key][c] + R * T * pylab.log(self.get_concentration(cids[c], c0, media))
            dG0_r[plot_key] = pylab.dot(S, dG0_f[plot_key])
            dG_r[plot_key] = pylab.dot(S, dG_f[plot_key])
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 2
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        profile_fig = pylab.figure()
        pylab.hold(True)
        data = pylab.zeros((Nr + 1, len(legend)))
        for i in range(len(legend)):
            for r in range(1, Nr + 1):
                data[r, i] = sum(dG_r[legend[i]][:r, 0])
        pylab.plot(data)
        legend(legend, loc="lower left")

        for i in range(len(rids)):
            pylab.text(i + 0.5, pylab.mean(data[i:(i + 2), 0]), rids[i], fontsize=6, horizontalalignment='center', backgroundcolor='white')
        
        pylab.xlabel("Reaction no.")
        pylab.ylabel("dG [kJ/mol]")
        self.html_writer.embed_matplotlib_figure(profile_fig, width=800, heigh=600)
        self.html_writer.write('</p>')
    
    def analyze_slack(self, key, field_map):
        self.html_writer.write('<p>\n')
        self.html_writer.write('<ul>\n')
        self.html_writer.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc) 
        c_mid = field_map.GetFloatField('C_MID', default_value=1e-3)
        (pH, I, T) = (default_pH, default_I, default_T)
        concentration_bounds = copy.deepcopy(self.kegg.cid2bounds)
        if ("CONDITIONS" in field_map):
            pH = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            I = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            T = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            self.html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                concentration_bounds[cid] = (conc, conc)
                self.html_writer.write(', [C%05d] = %g\n' % (cid, conc))
            self.html_writer.write('</li>\n')
        self.html_writer.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.kegg.write_reactions_to_html(self.html_writer, S, rids, fluxes, cids, show_cids=False)
        self.html_writer.write('</ul>\n')
        self.write_metabolic_graph(key, S, rids, cids)
        
        physiological_pC = field_map.GetFloatField('PHYSIO', default_value=4)
        (Nr, Nc) = S.shape

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = self.thermo.write_pseudoisomers_to_html(self.html_writer, self.kegg, cids)
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        pC = pylab.arange(0, 20, 0.1)
        B_vec = pylab.zeros(len(pC))
        #label_vec = [""] * len(pC)
        #limiting_reactions = set()
        for i in xrange(len(pC)):
            c_range = pC_to_range(pC[i], c_mid=c_mid)
            unused_dG_f, unused_concentrations, B = find_mtdf(S, dG0_f, c_range=c_range, bounds=bounds)
            B_vec[i] = B
            #curr_limiting_reactions = set(find(abs(dG_r - B) < 1e-9)).difference(limiting_reactions)
            #label_vec[i] = ", ".join(["%d" % rids[r] for r in curr_limiting_reactions]) # all RIDs of reactions that have dG_r = B
            #limiting_reactions |= curr_limiting_reactions

        try:
            unused_dG_f, unused_concentrations, pCr = find_pCr(S, dG0_f, c_mid=c_mid, bounds=bounds)
        except LinProgNoSolutionException:
            pCr = None
            
        try:
            c_range = pC_to_range(physiological_pC, c_mid=c_mid)
            unused_dG_f, unused_concentrations, B_physiological = find_mtdf(S, dG0_f, c_range=c_range, bounds=bounds) 
        except LinProgNoSolutionException:
            B_physiological = None

        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 5
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        slack_fig = pylab.figure()
        pylab.plot(pC, B_vec, 'b')
        #for i in xrange(len(pC)):
        #    text(pC[i], B_vec[i], label_vec[i], fontsize=6, horizontalalignment='left', backgroundcolor='white')
        pylab.xlabel('pC')
        pylab.ylabel('slack [kJ/mol]')
        (ymin, _) = pylab.ylim()
        (xmin, _) = pylab.xlim()
#        broken_barh([(xmin, pCr), (pCr, xmax)], (ymin, 0), facecolors=('yellow', 'green'), alpha=0.3)
        pylab.axhspan(ymin, 0, facecolor='b', alpha=0.15)
        title = 'C_mid = %g' % c_mid
        if (pCr != None and pCr < pC.max()):
            title += ', pCr = %.1f' % pCr
            pylab.plot([pCr, pCr], [ymin, 0], 'k--')
            pylab.text(pCr, 0, 'pCr = %.1f' % pCr, fontsize=8)
            if (pCr < physiological_pC):
                pylab.axvspan(pCr, physiological_pC, facecolor='g', alpha=0.3)
        if (B_physiological != None and physiological_pC < pC.max()):
            title += ', slack = %.1f [kJ/mol]' % B_physiological
            pylab.plot([xmin, physiological_pC], [B_physiological, B_physiological], 'k--')
            pylab.text(physiological_pC, B_physiological, 'B=%.1f' % B_physiological, fontsize=8)
        
        pylab.title(title)
        pylab.ylim(ymin=ymin)
        self.html_writer.embed_matplotlib_figure(slack_fig, width=800, height=600)

        # write a table of the compounds and their dG0_f
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f' [kJ/mol]"))
        for c in range(Nc):
            compound = self.kegg.cid2compound(cids[c])
            cid_str = '<a href="%s">C%05d</a>' % (compound.get_link(), compound.cid)
            self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%.1f</td>\n' % (cid_str, compound.name, dG0_f[c, 0]))
        self.html_writer.write('</table><br>\n')
        
        # write a table of the reactions and their dG0_r
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG RID", "Reaction", "flux"))
        for r in range(Nr):
            rid_str = '<a href="http://www.genome.jp/dbget-bin/www_bget?rn:R%05d">R%05d</a>' % (rids[r], rids[r])
            spr = {}
            for c in pylab.find(S[r, :]):
                spr[cids[c]] = S[r, c]
            reaction_str = self.kegg.sparse_to_hypertext(spr)
            self.html_writer.write('<tr><td>%s</td><td>%s</td><td>%g</td>\n' % (rid_str, reaction_str, fluxes[r]))
        self.html_writer.write('</table><br>\n')
        
        self.html_writer.write('</p>\n')

    def analyze_margin(self, key, field_map):
        self.html_writer.write('<p>\n')
        (S, rids, fluxes, cids) = self.get_reactions(key, field_map)
        self.html_writer.write('<li>Conditions:</br><ol>\n')
                    
        # The method for how we are going to calculate the dG0
        method = field_map.GetStringField('METHOD', default_value='MILO')

        if method == "MILO":
            thermodynamics = self.thermo
        elif method == "HATZI":
            thermodynamics = self.thermo.hatzi
        else:
            raise Exception("Unknown dG evaluation method: " + method)
        
        if "CONDITIONS" in field_map:
            thermodynamics.pH = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            thermodynamics.I = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            thermodynamics.T = ThermodynamicAnalysis.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            self.html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (thermodynamics.pH, thermodynamics.I, thermodynamics.T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                thermodynamics.bounds[cid] = (conc, conc)
                self.html_writer.write(', [C%05d] = %g\n' % (cid, conc))
            self.html_writer.write('</li>\n')
        self.html_writer.write('</ol></li>')
        thermodynamics.c_range = field_map.GetVFloatField('C_RANGE', default_value=[1e-6, 1e-2])
        thermodynamics.c_mid = field_map.GetFloatField('C_MID', default_value=1e-3)
        
        thermodynamic_pathway_analysis(S, rids, fluxes, cids, thermodynamics, self.html_writer)
        kegg_utils.write_module_to_html(self.html_writer, S, rids, fluxes, cids)

    def analyze_contour(self, key, field_map):
        pH_list = field_map.GetVFloatField("PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I_list = field_map.GetVFloatField("I", [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        T = field_map.GetFloatField("T", default_T)
        most_abundant = field_map.GetBoolField("ABUNDANT", False)
        formula = field_map.GetStringField("REACTION")
        c0 = field_map.GetFloatField("C0", 1.0)
        media = field_map.GetStringField("MEDIA", "None")
        if (media == "None"):
            media = None
            
        sparse_reaction = self.kegg.formula_to_sparse(formula)
        dG_r = self.estimate_dG_reaction(sparse_reaction, pH_list, I_list, T, c0, media, most_abundant)
        contour_fig = pylab.figure()
        
        pH_meshlist, I_meshlist = pylab.meshgrid(pH_list, I_list)
        CS = pylab.contour(pH_meshlist.T, I_meshlist.T, dG_r)       
        pylab.clabel(CS, inline=1, fontsize=10)
        pylab.xlabel("pH")
        pylab.ylabel("Ionic Strength")
        self.html_writer.embed_matplotlib_figure(contour_fig, width=800, height=600)
        self.html_writer.write('<br>\n' + self.kegg.sparse_to_hypertext(sparse_reaction) + '<br>\n')
        self.html_writer.write('</p>')

    def analyze_protonation(self, key, field_map):
        pH_list = field_map.GetVFloatField("PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I = field_map.GetFloatField("I", default_I)
        T = field_map.GetFloatField("T", default_T)
        cid = field_map.GetStringField("COMPOUND")
        cid = int(cid[1:])
        pmatrix = self.thermo.cid2PseudoisomerMap(cid).ToMatrix()
        data = pylab.zeros((len(self.cid), len(pH_list)))
        for j in range(len(pH_list)):
            pH = pH_list[j]
            dG0_array = pylab.matrix([-transform(dG0, nH, z, pH, I, T) / (R * T) \
                                      for (nH, z, dG0) in self.cid])
            dG0_array = dG0_array - max(dG0_array)
            p_array = pylab.exp(dG0_array)
            p_array = p_array / sum(p_array)
            data[:, j] = p_array    
        
        protonation_fig = pylab.figure()
        pylab.plot(pH_list, data.T)
        prop = matplotlib.font_manager.FontProperties(size=10)
        name = self.kegg.cid2name(cid)
        pylab.legend(['%s [%d]' % (name, z) for (nH, z, dG0) in pmatrix], prop=prop)
        pylab.xlabel("pH")
        pylab.ylabel("Pseudoisomer proportion")
        self.html_writer.embed_matplotlib_figure(protonation_fig, width=800, height=600)
        self.html_writer.write('<table border="1">\n')
        self.html_writer.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for (nH, z, dG0) in pmatrix:
            self.html_writer.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        self.html_writer.write('</table>')
        self.html_writer.write('</p>')

if __name__ == "__main__":
    db = SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter('../res/thermodynamic_pathway_analysis.html')
    G = GroupContribution(db, html_writer=html_writer)
    G.init()
    G.read_compound_abundance('../data/thermodynamics/compound_abundance.csv')
    thermo_analyze = ThermodynamicAnalysis(db, html_writer, thermodynamics=G)
    thermo_analyze.analyze_pathway("../data/thermodynamics/pathways.txt")
    
