import re, logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rms_flat
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException,\
    KeggReactionNotBalancedException
from pygibbs.thermodynamics import default_T, MissingCompoundFormationEnergy,\
    PsuedoisomerTableThermodynamics
from toolbox.util import _mkdir, calc_rmse
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_constants import R, symbol_dr_G0_prime
from pygibbs.kegg_reaction import Reaction
from toolbox.database import SqliteDatabase
import csv
import pydot
from toolbox.plotting import binned_plot
from pygibbs.thermodynamic_errors import MissingReactionEnergy
from collections import defaultdict

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
                "%s cannot use this NIST reaction because it is missing information about K', pH or I"
                % row_dict['reference_id'])
        
        self.T = NistRowData.none_float(row_dict['T']) or default_T
        
        # if there is no information about Ionic strength of pMg, assume by
        # default that the concentration of ions is 0 (note that pMg = 14 
        # is effectively [Mg] = 0).
        self.I = NistRowData.none_float(row_dict['I']) or 0.0
        self.pMg = NistRowData.none_float(row_dict['pMg']) or 14.0 
        self.dG0_r = -R*self.T*np.log(self.K_tag)
        self.evaluation = row_dict['evaluation']
        self.url = row_dict['url']
        self.ref_id = row_dict['reference_id']
        self.ec = row_dict['ec']
        kegg_reaction = row_dict['kegg_reaction']
        if not kegg_reaction:
            raise NistMissingCrucialDataException(
                "%s: cannot use this NIST reaction because it couldn't be mapped to KEGG IDs"
                % row_dict['reference_id'])
        try:
            self.reaction = NistRowData.ParseReactionFormula(name, kegg_reaction) 
        except KeggParseException as e:
            raise NistMissingCrucialDataException(
                "%s: cannot use reaction \"%s\", because: %s" %
                (row_dict['reference_id'], kegg_reaction, str(e)))

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
        other.ec = self.ec
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
        """ parse a two-sided formula such as: 2 C00001 = C00002 + C00003 
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
                dG0_f_tag = thermodynamics.GetTransformedFormationEnergies(cid, 
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
        self.pH_range = None
        self.override_I = None
        self.override_pMg = None
        self.override_T = None
        self.FromDatabase()
        self.BalanceReactions()

    def FromDatabase(self):
        self.data = []
        self.cid2count = {}
        logging.info('Reading NIST reaction data from database ...')
        for i, row_dict in enumerate(self.db.DictReader('nist_equilibrium')):
            nist_row_data = NistRowData()
            try:
                nist_row_data.ReadFromDatabase('nist%05d' % i, row_dict)
                self.data.append(nist_row_data)
                for cid in nist_row_data.GetAllCids():
                    self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1
            except NistMissingCrucialDataException as e:
                logging.debug(str(e))
        logging.info('Total of %d rows read from the NIST database' % len(self.data))
        
    def BalanceReactions(self, balance_water=True):
        for row in self.data:
            try:
                row.reaction.Balance(balance_water)
            except KeggReactionNotBalancedException as e:
                raise Exception(str(e) + '\n' + str(row.reaction) + '\n' + row.url)

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
        fig = plt.figure()
        plt.title("Temperature histogram")
        plt.hist(T_list, np.arange(int(min(T_list)), int(max(T_list)+1), 2.5))
        plt.xlabel("Temperature (C)")
        plt.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_T')

        fig = plt.figure()
        plt.hist(pMg_list, np.arange(0, 10.1, 0.1))
        plt.title("pMg histogram")
        plt.xlabel("pMg")
        plt.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_pMg')

        fig = plt.figure()
        plt.hist(pH_list, np.arange(4, 11, 0.1))
        plt.title("pH histogram")
        plt.xlabel("pH")
        plt.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_pH')

        fig = plt.figure()
        plt.hist(I_list, np.arange(0, 1, 0.025))
        plt.title("Ionic Strength histogram")
        plt.xlabel("Ionic Strength [M]")
        plt.ylabel("No. of measurements")
        html_writer.embed_matplotlib_figure(fig, width=320, height=240, name='hist_I')

        # histogram of publication years
        fig = plt.figure()
        plt.hist(year_list, np.arange(1930, 2010, 5))
        plt.title("Year of publication histogram")
        plt.xlabel("Year of publication")
        plt.ylabel("No. of measurements")
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
        hist_a = np.zeros(N)
        hist_b = np.zeros(N)
        for cid, cnt in self.cid2count.iteritems():
            if cnt >= N:
                cnt = N-1
            if cid in alberty_cids:
                hist_a[cnt] += 1
            else:
                hist_b[cnt] += 1
        hist_a[0] = len(alberty_cids.difference(self.cid2count.keys()))
        
        fig = plt.figure()
        plt.rc('font', size=10)
        plt.hold(True)
        p1 = plt.bar(range(N), hist_a, color='b')
        p2 = plt.bar(range(N), hist_b, color='r', bottom=hist_a[0:N])
        plt.text(N-1, hist_a[N-1]+hist_b[N-1], '> %d' % (N-1), fontsize=10, horizontalalignment='right', verticalalignment='baseline')
        plt.title("Overlap with Alberty's database")
        plt.xlabel("N reactions")
        plt.ylabel("no. of compounds measured in N reactions")
        plt.legend((p1[0], p2[0]), ("Exist in Alberty's database", "New compounds"))

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

    def verify_formation(self, html_writer, thermodynamics, name=None):
        cid2errors = defaultdict(list)
        cid2refs = defaultdict(set)
        reaction2errors = defaultdict(list)
        reaction2refs = defaultdict(set)
        for row_data in self.SelectRowsFromNist():
            dG0_est = row_data.PredictReactionEnergy(thermodynamics)
            if np.isnan(dG0_est):
                continue
            err = row_data.dG0_r - dG0_est
            for cid in row_data.GetAllCids():
                cid2errors[cid].append(err)
                cid2refs[cid].add((row_data.ref_id, row_data.url))
            reaction2errors[row_data.reaction].append(err)
            reaction2refs[row_data.reaction].add((row_data.ref_id, row_data.url))
        
        rowdicts = []
        for cid, err_list in cid2errors.iteritems():
            refs = cid2refs[cid]
            urls = ', '.join(['<a href="%s">%s</a>' % (url, ref_id)
                              for ref_id, url in refs])
            rowdict = {'cid':'C%05d' % cid,
                       'name':self.kegg.cid2name(cid),
                       'RMSE':rms_flat(err_list),
                       'E[err]':np.mean(err_list),
                       '#err':len(err_list),
                       'std[err]':np.std(err_list),
                       'URLs':urls}
            rowdicts.append(rowdict)
        
        rowdicts.sort(key=lambda x:x['RMSE'], reverse=True)
        html_writer.write_table(rowdicts, ['#', 'cid', 'name', 'RMSE',
                                           '#err', 'E[err]', 'std[err]', 'URLs'], decimal=1)
        
        rowdicts = []
        for reaction, err_list in reaction2errors.iteritems():
            refs = reaction2refs[reaction]
            urls = ', '.join(['<a href="%s">%s</a>' % (url, ref_id)
                              for ref_id, url in refs])
            rowdict = {'reaction':reaction.to_hypertext(show_cids=False),
                       'RMSE':rms_flat(err_list),
                       'E[err]':np.mean(err_list),
                       '#err':len(err_list),
                       'std[err]':np.std(err_list),
                       'URLs':urls}
            rowdicts.append(rowdict)
        
        rowdicts.sort(key=lambda x:x['RMSE'], reverse=True)
        html_writer.write_table(rowdicts, ['#', 'reaction', 'RMSE',
                                           '#err', 'E[err]', 'std[err]', 'URLs'], decimal=1)
        
        
    def verify_results(self, html_writer, thermodynamics, name=None):
        """Calculate all the dG0_r for the reaction from NIST and compare to
           the measured data.
        
        Write results to HTML.
        
        Args:
            thermodynamics: a Thermodynamics object that provides dG estimates.
            ignore_I: whether or not to ignore the ionic strength in NIST.
        """
        
        dG0_obs_vec = []
        dG0_est_vec = []
       
        # A mapping from each evaluation method (NIST calls separates them to
        # A, B, C and D) to the results of the relevant measurements
        evaluation_map = {}
        rowdicts = []
        finite_rowdicts = []
        
        eval_to_label = {'A':'high quality', 'B':'low quality', 'C':'low quality', 'D':'low quality'}
        
        for row_data in self.SelectRowsFromNist():
            rowdict = {}
            label = eval_to_label[row_data.evaluation]
            if label not in evaluation_map:
                evaluation_map[label] = ([], [])
            rowdict[symbol_dr_G0_prime + ' (obs)'] = np.round(row_data.dG0_r, 1)
            rowdict['_reaction'] = row_data.reaction
            rowdict['reaction'] = row_data.reaction.to_hypertext(show_cids=False)
            if row_data.reaction.rid is not None:
                rowdict['rid'] = '<a href="%s">R%05d</a>' % (row_data.reaction.get_link(), row_data.reaction.rid)
            else:
                rowdict['rid'] = ''
            rowdict['pH'] = row_data.pH
            rowdict['pMg'] = row_data.pMg
            rowdict['I'] = row_data.I
            rowdict['T'] = row_data.T
            rowdict['eval.'] = row_data.evaluation
            rowdict['url'] = '<a href="%s">%s</a>' % (row_data.url, row_data.ref_id)

            dG0_est = row_data.PredictReactionEnergy(thermodynamics)
            if np.isfinite(dG0_est):
                dG0_obs_vec.append(row_data.dG0_r)
                dG0_est_vec.append(dG0_est)
                evaluation_map[label][0].append(row_data.dG0_r)
                evaluation_map[label][1].append(dG0_est)
                rowdict[symbol_dr_G0_prime + ' (est)'] = np.round(dG0_est, 1)
                rowdict['residual'] = np.round(row_data.dG0_r - dG0_est, 3)
                rowdict['|error|'] = abs(rowdict['residual'])
                rowdict['sort_key'] = -rowdict['|error|']
                finite_rowdicts.append(rowdict)
            else:
                rowdict['sort_key'] = 1
            
            rowdicts.append(rowdict)
        
        rowdicts.sort(key=lambda x:x['sort_key'])
        
        if not dG0_obs_vec:
            return 0, 0

        unique_reaction_dict = defaultdict(list)
        for rowdict in finite_rowdicts:
            unique_reaction_dict[rowdict['_reaction']].append(rowdict['|error|'])
        unique_rmse_list = [rms_flat(error_list)
                            for error_list in unique_reaction_dict.values()]
        unique_rmse = rms_flat(unique_rmse_list)
        
        rmse = calc_rmse(dG0_obs_vec, dG0_est_vec)

        # plot the profile graph
        plt.rcParams['text.usetex'] = False
        plt.rcParams['legend.fontsize'] = 10
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.size'] = 12
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['lines.markersize'] = 3
        
        fig1 = plt.figure(figsize=(6,6), dpi=90)
        plt.hold(True)
        
        colors = ['purple', 'orange']
        for i, label in enumerate(sorted(evaluation_map.keys())):
            measured, predicted = evaluation_map[label]
            plt.plot(measured, predicted, marker='.', linestyle='None', 
                       markerfacecolor=colors[i], markeredgecolor=colors[i], 
                       markersize=5, label=label, figure=fig1)
        
        plt.legend(loc='lower right')
        
        plt.text(-50, 40, r'RMSE = %.1f [kJ/mol]' % (unique_rmse), fontsize=14,
                 figure=fig1)
        plt.xlabel(r'observed $\Delta_r G^{\'\circ}$ [kJ/mol]', fontsize=14, figure=fig1)
        plt.ylabel(r'estimated $\Delta_r G^{\'\circ}$ [kJ/mol]', fontsize=14, figure=fig1)
        #min_x = min(dG0_obs_vec)
        #max_x = max(dG0_obs_vec)
        plt.plot([-60, 60], [-60, 60], 'k--', figure=fig1)
        plt.axis([-60, 60, -60, 60])
        if name:
            html_writer.embed_matplotlib_figure(fig1, name=name+"_eval")
        else:
            html_writer.embed_matplotlib_figure(fig1)
        
        fig2 = plt.figure(figsize=(6,6), dpi=90)
        binned_plot(x=[rowdict['pH'] for rowdict in finite_rowdicts],
                    y=[rowdict['|error|'] for rowdict in finite_rowdicts],
                    bins=[5,6,7,8,9],
                    y_type='rmse',
                    figure=fig2)
        plt.xlim((4, 11))
        plt.ylim((0, 12))
        plt.title(r'effect of pH', fontsize=14, figure=fig2)
        plt.xlabel('pH', fontsize=14, figure=fig2)
        plt.ylabel(r'RMSE ($\Delta_r G^{\'\circ}$) [kJ/mol]', 
                   fontsize=14, figure=fig2)
        if name:
            html_writer.embed_matplotlib_figure(fig2, name=name+"_pH")
        else:
            html_writer.embed_matplotlib_figure(fig2)
        
        fig3 = plt.figure(figsize=(6,6), dpi=90)
        plt.hist([rowdict['residual'] for rowdict in finite_rowdicts],
                 bins=np.arange(-50, 50, 0.5))
        plt.title(r'RMSE = %.1f [kJ/mol]' % rmse, fontsize=14, figure=fig3)
        plt.xlabel(r'residual $\Delta_r G^{\'\circ}$ [kJ/mol]',
                   fontsize=14, figure=fig3)
        plt.ylabel(r'no. of measurements', fontsize=14, figure=fig3)
        if name:
            html_writer.embed_matplotlib_figure(fig3, name=name+"_hist")
        else:
            html_writer.embed_matplotlib_figure(fig3)

        table_headers = ["#", "|error|",
                         symbol_dr_G0_prime + " (obs)",
                         symbol_dr_G0_prime + " (est)",
                         "reaction", "rid", "pH", "pMg", "I", "T",
                         "eval.", "url"]
        html_writer.write_table(rowdicts, table_headers, decimal=1)
        
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
            except MissingReactionEnergy as e:
                logging.debug("the reaction in (%s) cannot be estimated: %s" % (row_data.ref_id, str(e)))
                continue
                
            total_list.append([row_data.dG0_r, dG0_pred1, dG0_pred2, 
                               row_data.reaction, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               row_data.url])
        
        if not total_list:
            return 0, 0
        
        # plot the profile graph
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.size'] = 8
        plt.rcParams['lines.linewidth'] = 2
        plt.rcParams['lines.markersize'] = 2
        plt.rcParams['figure.dpi'] = 100
        
        data_mat = np.array(total_list)
        fig1 = plt.figure(figsize=(4,4))
        plt.hold(True)
        error1 = data_mat[:,0]-data_mat[:,1]
        error2 = data_mat[:,0]-data_mat[:,2]
        
        max_err = max(error1.max(), error2.max())
        min_err = min(error1.min(), error2.min())
        plt.plot([min_err, max_err], [min_err, max_err], 'k--', figure=fig1)
        plt.plot(error1, error2, '.', figure=fig1)
        plt.title("Error Comparison per Reaction (in kJ/mol)")
        plt.xlabel(thermo1.name, figure=fig1)
        plt.ylabel(thermo2.name, figure=fig1)
        html_writer.embed_matplotlib_figure(fig1, name=name+"_corr")
        
        fig2 = plt.figure(figsize=(7,3))
        for i, thermo in enumerate([thermo1, thermo2]):
            fig2.add_subplot(1,2,i+1)
            plt.plot(data_mat[:,0], data_mat[:,i+1], 'b.')
            rmse = calc_rmse(data_mat[:,0], data_mat[:,i+1])
            plt.text(-50, 40, r'RMSE = %.1f [kJ/mol]' % (rmse))
            plt.xlabel(r'observed $\Delta G_r^\circ$ from NIST [kJ/mol]')
            plt.ylabel(r'estimated $\Delta G_r^\circ$ using %s [kJ/mol]' % thermo.name)
            plt.plot([-60, 60], [-60, 60], 'k--')
            plt.axis([-60, 60, -60, 60])
        
        html_writer.embed_matplotlib_figure(fig2, name=name+"_eval")

        table_headers = ["dG'0 (obs)", "dG'0 (%s)" % thermo1.name, 
                         "dG'0 (%s)" % thermo2.name, "reaction", "rid", "pH", 
                         "pMg", "I", "T", "eval.", "url"]
        dict_list = []
        for row in sorted(total_list, key=lambda(x):abs(x[1]-x[2]), reverse=True):
            d = {}
            d["dG'0 (obs)"] = '%.1f' % row[0]
            d["dG'0 (%s)" % thermo1.name] = '%.1f' % row[1]
            d["dG'0 (%s)" % thermo2.name] = '%.1f' % row[2]
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
            
    def SelectRowsFromNist(self, reaction=None, check_reverse=True, 
                           T_range=None, pH_range=None):
        T_range = T_range or self.T_range
        pH_range = pH_range or self.pH_range
        rows = []
        checklist = []
        if reaction:
            checklist.append(reaction)
            if check_reverse:
                checklist.append(reaction.reverse())
        for nist_row_data in self.data:
            if T_range and not (T_range[0] < nist_row_data.T < T_range[1]):
                continue 
            if pH_range and not (pH_range[0] < nist_row_data.pH < pH_range[1]):
                continue 
            if checklist and nist_row_data.reaction not in checklist:
                continue
            if self.override_pMg or self.override_I or self.override_T:
                nist_row_copy = nist_row_data.Clone()
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
    fp = open('../res/nist_kegg_ids.txt', 'w')
    for cid in nist.GetAllCids():
        fp.write("C%05d\n" % cid)
    fp.close()
    nist.AnalyzeStats(html_writer)
    nist.AnalyzeConnectivity(html_writer)
    html_writer.close()
