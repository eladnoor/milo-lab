import pylab, re, logging
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException
from pygibbs.thermodynamics import default_T, MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from toolbox.util import _mkdir, calc_rmse
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_constants import R
from pygibbs.kegg_reaction import Reaction
from toolbox.database import SqliteDatabase
import copy
import csv
import pydot
from toolbox.plotting import binned_plot

class NistMissingCrucialDataException(Exception):
    pass


class NistReactionBalanceException(Exception):
    pass


class NistRowData:
    def __init__(self):
        pass
    
    @staticmethod
    def none_float(x):
        if x:
            return float(x)
        else:
            return None
    
    def ReadFromDatabase(self, name, row_dict):
        self.K_tag = NistRowData.none_float(row_dict['K_tag'])
        self.pH = NistRowData.none_float(row_dict['pH'])
        if not (self.K_tag and self.pH): # missing Keq or pH makes the data unusable
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it is missing information about K', pH or I")
        
        self.T = NistRowData.none_float(row_dict['T']) or default_T
        
        # if there is no imformation about Ionic strength of pMg, assume by
        # default that the concentration of ions is 0 (note that pMg = 14 
        # is effectively [Mg] = 0).
        self.I = NistRowData.none_float(row_dict['I']) or 0.0
        self.pMg = NistRowData.none_float(row_dict['pMg']) or 14.0 
        self.dG0_r = -R*self.T*pylab.log(self.K_tag)
        self.evaluation = row_dict['evaluation']
        self.url = row_dict['url']
        self.ref_id = row_dict['reference_id']
        kegg_reaction = row_dict['kegg_reaction']
        if not kegg_reaction:
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it couldn't be mapped to KEGG IDs")
        try:
            self.reaction = NistRowData.ParseReactionFormula(name, kegg_reaction) 
        except KeggParseException as e:
            raise NistMissingCrucialDataException("cannot use reaction \"%s\", because: %s" % (kegg_reaction, str(e)))

    def Clone(self):
        other = NistRowData()
        other.K_tag = self.K_tag
        other.pH = self.pH
        other.I = self.I
        other.T = self.T
        other.pMg = self.pMg
        other.dG0_r = self.dG0_r
        other.evaluation = self.evaluation
        other.url = self.url
        other.ref_id = self.ref_id
        other.reaction = self.reaction.clone()
        return other

    def GetYear(self):
        try:
            year = int(self.ref_id[0:2])
            if year < 11:
                return year + 2000
            else:
                return year + 1900
        except ValueError:
            return None
    
    @staticmethod
    def ParseReactionFormula(name, formula):
        """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
            return the set of substrates, products and the direction of the reaction
        """
        try:
            left, right = formula.split(' = ', 1)
        except ValueError:
            raise KeggParseException("There should be exactly one '=' sign")
        sparse_reaction = {}
        for cid, amount in NistRowData.ParseReactionFormulaSide(left).iteritems():
            sparse_reaction[cid] = -amount
        for cid, amount in NistRowData.ParseReactionFormulaSide(right).iteritems():
            if (cid in sparse_reaction):
                raise KeggParseException("C%05d appears on both sides of this formula" % cid)
            sparse_reaction[cid] = amount
        
        reaction = Reaction([name], sparse_reaction, None, '=>')
        
        kegg = Kegg.getInstance()
        rid = kegg.reaction2rid(reaction) or kegg.reaction2rid(reaction.reverse())
        reaction.rid = rid
        return reaction
    
    @staticmethod
    def ParseReactionFormulaSide(s):
        """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            return the set of CIDs, ignore stoichiometry
        """
        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if (len(tokens) == 1):
                amount = 1
                cid = member
            else:
                amount = tokens[0]
                cid = tokens[1]
            if cid[0] != 'C':
                raise KeggParseException("Compound ID does not start with a C: " + cid)
            try:
                cid = int(cid[1:])
            except ValueError:
                raise KeggParseException("Compound ID is not an integer number: " + cid)
            compound_bag[cid] = compound_bag.get(cid, 0) + float(amount)
        
        return compound_bag
    
    def GetAllCids(self):
        return set(self.reaction.get_cids())
    
    def PredictReactionEnergy(self, thermodynamics):
        return self.reaction.PredictReactionEnergy(thermodynamics, pH=self.pH,
            pMg=self.pMg, I=self.I, T=self.T)
            
    def NormalizeCompounds(self, thermodynamics):
        """
            Subtracts the formation energies of known compounds (i.e. those
            that appear in 'thermodynamics'), and removes them from the reaction.
            This is useful for algorithms such as linear regression where one
            would want to fix the values for some of the compounds and estimate
            only those of the unknown ones.
        """
        cids_to_remove = set()
        for cid, coeff in self.reaction.sparse.iteritems():
            try:
                dG0_f_tag = thermodynamics.cid2dG0_tag(cid, 
                    pH=self.pH, pMg=self.pMg, I=self.I, T=self.T)
                cids_to_remove.add(cid)
                self.dG0_r -= coeff * dG0_f_tag
            except MissingCompoundFormationEnergy:
                continue
        
        for cid in cids_to_remove:
            del self.reaction.sparse[cid]
    
class Nist(object):
    def __init__(self, T_range=(298, 314)):
        self.db = SqliteDatabase('../data/public_data.sqlite')
        self.kegg = Kegg.getInstance()
        self.T_range = T_range
        self.override_I = None
        self.override_pMg = None
        self.override_T = None
        self.FromDatabase()
        self.BalanceReactions()

    def FromDatabase(self):
        self.data = []
        self.cid2count = {}
        logging.info('Reading NIST reaction data from database')
        for i, row_dict in enumerate(self.db.DictReader('nist_equilibrium')):
            nist_row_data = NistRowData()
            try:
                nist_row_data.ReadFromDatabase('nist%05d' % i, row_dict)
                self.data.append(nist_row_data)
                for cid in nist_row_data.GetAllCids():
                    self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1
            except NistMissingCrucialDataException:
                continue
        
    def BalanceReactions(self, balance_water=True):
        for row in self.data:
            row.reaction.Balance(balance_water)

    def GetAllCids(self):
        return sorted(self.cid2count.keys())
    
    def AnalyzeStats(self, html_writer):
        """
            Produces a set of plots that show some statistics about the NIST database
        """
        logging.info('Calculating statistics for NIST database (%d rows)' % len(self.data))
        
        if not self.data:
            raise Exception("The database has no rows in it")
        
        T_list = []
        I_list = []
        pH_list = []
        pMg_list = []
        year_list = []
        for nist_row_data in self.data:
            pH_list.append(nist_row_data.pH)
            T_list.append(nist_row_data.T - 273.15)

            if nist_row_data.I:
                I_list.append(nist_row_data.I)
            if nist_row_data.pMg:
                pMg_list.append(nist_row_data.pMg)
            
            year = nist_row_data.GetYear()
            if year:
                year_list.append(year)
        
        html_writer.write("<p><h2>NIST database statistics</h2>\n")
        fig = pylab.figure()
        pylab.title("Temperature histogram")
        pylab.hist(T_list, pylab.arange(int(min(T_list)), int(max(T_list)+1), 2.5))
        pylab.xlabel("Temperature (C)")
        pylab.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_T')

        fig = pylab.figure()
        pylab.hist(pMg_list, pylab.arange(0, 10.1, 0.1))
        pylab.title("pMg histogram")
        pylab.xlabel("pMg")
        pylab.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_pMg')

        fig = pylab.figure()
        pylab.hist(pH_list, pylab.arange(4, 11, 0.1))
        pylab.title("pH histogram")
        pylab.xlabel("pH")
        pylab.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_pH')

        fig = pylab.figure()
        pylab.hist(I_list, pylab.arange(0, 1, 0.025))
        pylab.title("Ionic Strength histogram")
        pylab.xlabel("Ionic Strength [M]")
        pylab.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_I')

        # histogram of publication years
        fig = pylab.figure()
        pylab.hist(year_list, pylab.arange(1930, 2010, 5))
        pylab.title("Year of publication histogram")
        pylab.xlabel("Year of publication")
        pylab.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_year')

        db_public = SqliteDatabase('../data/public_data.sqlite')
        alberty = PsuedoisomerTableThermodynamics.FromDatabase(db_public, 'alberty_pseudoisomers')
        alberty_cids = set(alberty.get_all_cids())
        nist_cids = set(self.GetAllCids())
        
        count_list = ["Alberty #compounds = %d" % len(alberty_cids),
                      "NIST #compounds = %d" % len(nist_cids),
                      "intersection #compounds = %d" % len(alberty_cids.intersection(nist_cids))]
        html_writer.write_ul(count_list)
        
        N = 60 # cutoff for the number of counts in the histogram
        hist_a = pylab.zeros(N)
        hist_b = pylab.zeros(N)
        for cid, cnt in self.cid2count.iteritems():
            if cnt >= N:
                cnt = N-1
            if cid in alberty_cids:
                hist_a[cnt] += 1
            else:
                hist_b[cnt] += 1
        hist_a[0] = len(alberty_cids.difference(self.cid2count.keys()))
        
        fig = pylab.figure()
        pylab.rc('font', size=10)
        pylab.hold(True)
        p1 = pylab.bar(range(N), hist_a, color='b')
        p2 = pylab.bar(range(N), hist_b, color='r', bottom=hist_a[0:N])
        pylab.text(N-1, hist_a[N-1]+hist_b[N-1], '> %d' % (N-1), fontsize=10, horizontalalignment='right', verticalalignment='baseline')
        pylab.title("Overlap with Alberty's database")
        pylab.xlabel("N reactions")
        pylab.ylabel("no. of compounds measured in N reactions")
        pylab.legend((p1[0], p2[0]), ("Exist in Alberty's database", "New compounds"))

        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='connectivity')

    def AnalyzeConnectivity(self, html_writer):
        
        def cid2name(cid, KEGG):
            return "\"" + KEGG.cid2name(cid) + "\""
        
        def load_cid_set(train_csv_fname):
            """
                Read the training data from a CSV file
            """
            cid_set = set()
            for row in csv.DictReader(open(train_csv_fname)):
                #(smiles, cid, compoud_name, dG0, dH0, z, nH, Mg, use_for, ref, remark) = row
                if (row['use for'] in ['skip']):
                    continue
                cid = int(row['cid'])
                if cid > 0:
                    cid_set.add(cid)
            return cid_set
                
        known_cids = load_cid_set('../data/thermodynamics/dG0_seed.csv')
        one_step_cids = set()
        coupled_cids = set()
        Gdot = pydot.Dot()
        
        for nist_row_data in nist.data:
            unknown_cids = list(nist_row_data.GetAllCids().difference(known_cids))
            if len(unknown_cids) == 1:
                one_step_cids.add(unknown_cids[0])
            elif len(unknown_cids) == 2:
                coupled_cids.add((min(unknown_cids), max(unknown_cids)))
        
        for cid in one_step_cids:
            #Gdot.add_node(pydot.Node(cid2name(cid, KEGG), None))
            Gdot.add_node(pydot.Node("C%05d" % cid, None))
        
        for (cid1, cid2) in coupled_cids:
            Gdot.add_edge(pydot.Edge("C%05d" % cid1, "C%05d" % cid2, None))
        
        html_writer.write("<p><h2>Connectivity</h2>\n")
        html_writer.embed_dot_inline(Gdot, width=640, height=480)
        html_writer.write("</p>\n")
        #win = xdot.DotWindow()
        #win.connect('destroy', gtk.main_quit)
        #win.set_filter('dot')
        #util._mkdir('../res/nist')
        #dot_fname = '../res/nist/connectivity.dot'
        #Gdot.write(dot_fname, format='dot')
        #win.open_file(dot_fname)
        #gtk.main()

    def verify_results(self, html_writer, thermodynamics, name=None):
        """Calculate all the dG0_r for the reaction from NIST and compare to
           the measured data.
        
        Write results to HTML.
        
        Args:
            thermodynamics: a Thermodynamics object that provides dG estimates.
            ignore_I: whether or not to ignore the ionic strength in NIST.
        """
        
        known_cid_set = thermodynamics.get_all_cids()
        dG0_obs_vec = []
        dG0_est_vec = []
       
        # A mapping from each evaluation method (NIST calls separates them to
        # A, B, C and D) to the results of the relevant measurements
        evaluation_map = {}
        total_list = []
        
        eval_to_label = {'A':'high quality', 'B':'low quality', 'C':'low quality', 'D':'low quality'}
        
        for row_data in self.SelectRowsFromNist():
            row_cids = set(row_data.GetAllCids())
            unknown_cids = row_cids.difference(known_cid_set)
            if unknown_cids:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.ref_id)
                continue
            
            label = eval_to_label[row_data.evaluation]
            
            if label not in evaluation_map:
                evaluation_map[label] = ([], [])
            
            try:
                dG0_pred = row_data.PredictReactionEnergy(thermodynamics)
            except MissingCompoundFormationEnergy:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.origin)
                continue
                
            dG0_obs_vec.append(row_data.dG0_r)
            dG0_est_vec.append(dG0_pred)
            evaluation_map[label][0].append(row_data.dG0_r)
            evaluation_map[label][1].append(dG0_pred)
            error = abs(row_data.dG0_r - dG0_pred)

            total_list.append([error, row_data.dG0_r, dG0_pred, 
                               row_data.reaction, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               row_data.url])
        
        if not dG0_obs_vec:
            return 0, 0

        unique_reaction_dict = {}
        for error, _dG0_obs, _dG0_est, reaction, _pH, _pMg, _I, _T, _eval, _url in total_list: 
            unique_reaction_dict.setdefault(reaction, []).append((error)**2)
        unique_mean_sqr_error = [pylab.mean(sqr_error_list) for sqr_error_list in unique_reaction_dict.values()]
        unique_rmse = pylab.sqrt(pylab.mean(unique_mean_sqr_error))
        
        rmse = calc_rmse(dG0_obs_vec, dG0_est_vec)

        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 1
        pylab.rcParams['lines.markersize'] = 3
        pylab.rcParams['figure.figsize'] = [6.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        fig = pylab.figure()
        pylab.hold(True)
        
        colors = ['purple', 'orange']
        for i, label in enumerate(sorted(evaluation_map.keys())):
            measured, predicted = evaluation_map[label]
            pylab.plot(measured, predicted, marker='.', linestyle='None', 
                       markerfacecolor=colors[i], markeredgecolor=colors[i], 
                       markersize=5, label=label)
        
        pylab.legend(loc='lower right')
        
        pylab.text(-50, 40, r'RMSE = %.1f [kJ/mol]' % (unique_rmse), fontsize=14)
        pylab.xlabel(r'observed $\Delta G_r^\circ$ [kJ/mol]', fontsize=14)
        pylab.ylabel(r'estimated $\Delta G_r^\circ$ [kJ/mol]', fontsize=14)
        #min_x = min(dG0_obs_vec)
        #max_x = max(dG0_obs_vec)
        pylab.plot([-60, 60], [-60, 60], 'k--')
        pylab.axis([-60, 60, -60, 60])
        if name:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300, name=name+"_eval")
        else:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300)
        
        fig = pylab.figure()
        binned_plot(x=[row[4] for row in total_list], # pH
                    y=[row[0] for row in total_list],
                    bins=[6,8],
                    y_type='rmse',
                    figure=fig)
        pylab.xlim((4, 11))
        pylab.ylim((0, 12))
        pylab.title(r'effect of pH', fontsize=14, figure=fig)
        pylab.xlabel('pH', fontsize=14, figure=fig)
        pylab.ylabel(r'average $|\Delta_{obs} G^\circ - \Delta_{est} G^\circ|$ [kJ/mol]', fontsize=14, figure=fig)
        if name:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300, name=name+"_pH")
        else:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300)
        
        fig = pylab.figure()
        pylab.hist([(row[1] - row[2]) for row in total_list], bins=pylab.arange(-50, 50, 0.5))
        pylab.title(r'RMSE = %.1f [kJ/mol]' % rmse, fontsize=14)
        pylab.xlabel(r'$\Delta_{obs} G^\circ - \Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
        pylab.ylabel(r'no. of measurements', fontsize=14)
        if name:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300, name=name+"_hist")
        else:
            html_writer.embed_matplotlib_figure(fig, width=400, height=300)

        table_headers = ["|err|", "dG'0 (obs)", "dG'0 (est)", "reaction", "rid", "pH", "pMg", "I", "T", "eval.", "url"]
        dict_list = []
        for row in sorted(total_list, reverse=True):
            d = {}
            d['|err|'] = '%.1f' % row[0]
            d['dG\'0 (obs)'] = '%.1f' % row[1]
            d['dG\'0 (est)'] = '%.1f' % row[2]
            d['reaction'] = row[3].to_hypertext(show_cids=False)
            if row[3].rid is not None:
                d['rid'] = '<a href="%s">R%05d</a>' % (row[3].get_link(), row[3].rid)
            else:
                d['rid'] = ''
            d['pH'] = '%.1f' % row[4]
            d['pMg'] = '%.1f' % row[5]
            d['I'] = '%.2f' % row[6]
            d['T'] = '%.1f' % row[7]
            d['eval.'] = row[8]
            if row[9]:
                d['url'] = '<a href="%s">link</a>' % row[9]
            else:
                d['url'] = ''
            dict_list.append(d)
        html_writer.write_table(dict_list, table_headers)
        
        return len(dG0_obs_vec), unique_rmse
    
    def two_way_comparison(self, html_writer, thermo1, thermo2, name=None):
        """
            Compare the estimation errors of two different evaluation methods.
        
        Write results to HTML.
        
        Args:
            thermo1: a Thermodynamics object that provides dG estimates.
            thermo2: a Thermodynamics object that provides dG estimates.
        """
        
        total_list = []
        
        for row_data in self.SelectRowsFromNist():
            try:
                dG0_pred1 = row_data.PredictReactionEnergy(thermo1)
                dG0_pred2 = row_data.PredictReactionEnergy(thermo2)
            except MissingCompoundFormationEnergy:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.ref_id)
                continue
                
            total_list.append([row_data.dG0_r, dG0_pred1, dG0_pred2, 
                               row_data.reaction, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               row_data.url])
        
        if not total_list:
            return 0, 0
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 12
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 3
        pylab.rcParams['figure.figsize'] = [6.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        data_mat = pylab.array(total_list)
        fig = pylab.figure()
        pylab.hold(True)
        error1 = data_mat[:,0]-data_mat[:,1]
        error2 = data_mat[:,0]-data_mat[:,2]
        
        max_err = max(error1.max(), error2.max())
        min_err = min(error1.min(), error2.min())
        pylab.plot([min_err, max_err], [min_err, max_err], 'k--', figure=fig)
        pylab.plot(error1, error2, '.', figure=fig)
        pylab.title("Error Comparison per Reaction (in kJ/mol)")
        pylab.xlabel(thermo1.name, figure=fig)
        pylab.ylabel(thermo2.name, figure=fig)
        html_writer.embed_matplotlib_figure(fig, width=200, height=200, name=name+"_corr")

    def SelectRowsFromNist(self, reaction=None, check_reverse=True):
        rows = []
        checklist = []
        if reaction:
            checklist.append(reaction)
            if check_reverse:
                checklist.append(reaction.reverse())
        for nist_row_data in self.data:
            if self.T_range and not (self.T_range[0] < nist_row_data.T < self.T_range[1]):
                continue 
            if checklist and nist_row_data.reaction not in checklist:
                continue
            if self.override_pMg or self.override_I or self.override_T:
                nist_row_copy = copy.deepcopy(nist_row_data)
                if self.override_pMg:
                    nist_row_copy.pMg = self.override_pMg
                if self.override_I:
                    nist_row_copy.I = self.override_I
                if self.override_T:
                    nist_row_copy.T = self.override_T
                rows.append(nist_row_copy)
            else:
                rows.append(nist_row_data)
        return rows
    
    def GetUniqueReactionSet(self):
        return set([row.reaction for row in self.data])


if __name__ == '__main__':
    #logging.getLogger('').setLevel(logging.DEBUG)
    _mkdir("../res/nist")
    html_writer = HtmlWriter("../res/nist/statistics.html")
    nist = Nist()
    nist.AnalyzeStats(html_writer)
    nist.AnalyzeConnectivity(html_writer)
    html_writer.close()
