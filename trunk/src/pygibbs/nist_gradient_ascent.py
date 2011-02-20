import csv
from pylab import * #@UnusedWildImport
from toolbox.html_writer import HtmlWriter
from pygibbs.groups import GroupContribution
from pygibbs.kegg_errors import KeggParseException
from pygibbs.thermodynamic_constants import R
from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.alberty import Alberty
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from pygibbs import reverse_thermo
from pygibbs import pseudoisomer
from toolbox import util, database
import logging
from copy import deepcopy
from pygibbs.nist_regression import NistAnchors

class GradientAscent(Thermodynamics):
    def __init__(self, gc):
        Thermodynamics.__init__(self)
        self.gc = gc
        self.kegg = gc.kegg()
        self.data = []           # this will hold all the training data in one place, each row represents a measurement
        self.train_rowids = []  # the row indices in self.data to use for training
        self.test_rowids = []   # the row indices in self.data to use for testing
        self.cid2rowids = {} # holds a list of indices to self.train_data where the CID is participating
        self.cid2pmap_dict = {}       # this is the solution vector (what we are trying to solve)
        self.cache_error = {}    # for each rowid - holds the last calculated squared error (the difference between calc and measured)
        self.cache_cid = {}      # for each (rowid,cid) pair - holds the last calculated dG_f (multiplied by the stoichiometric coeff)

    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

    def load_dG0_data(self, train_csv_fname):
        """
            Read the training data from a CSV file, into the cid2pmap_dict.
            Returns a set containing the imported CIDs.
        """
        nist_anchors = NistAnchors(self.gc.db, self.html_writer)
        cids = set()
        for row_dict in csv.DictReader(open(train_csv_fname)):
            if row_dict.get('use_for', 'skip') == 'skip':
                continue
            
            cid = int(row_dict['cid'])
            dG0 = float(row_dict['dG0'])
            z = int(row_dict['z'])
            nH = int(row_dict['nH'])
            nMg = int(row_dict.get('nMg', 0))
            
            if cid:
                self.cid2pmap_dict.setdefault(cid, pseudoisomer.PseudoisomerMap())
                self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)
                cids.add(cid)
                self.cid2source_string = 'Gradient Ascent'
        return cids

    def load_nist_data(self, nist, thermodynamics, override_I=None, skip_missing_reactions=False, T_range=None):
        """
            skip_missing_reactions - if True, the method will not add reactions which have a compound
            with a no data (on dG0_f). In this case, the GradientAscent object is meant only to improve
            the dG0_f of the training data, but not to discover formation energies of new compounds.
        """
        self.cid2rowids = {}
        
        known_cids = set(thermodynamics.get_all_cids())
        
        # filter out the rows in nist.data that cannot be evaluated for some reason
        self.data = []
        for nist_row_data in nist.data:
            row_data = deepcopy(nist_row_data)
            if T_range and not (T_range[0] < row_data.T < T_range[1]):
                logging.warning('Temperature %f not within allowed range.', nist_row_data.T)
                continue # the temperature is outside the allowed range
            
            if override_I:
                row_data.I = override_I
                
            if row_data.evaluation != 'A':
                row_data.evaluation = 'B - D'
                
            #evaluation = row[2]
            reaction_cids = set(row_data.GetAllCids())
            unknown_cids = reaction_cids.difference(known_cids)
            
            if unknown_cids:
                if skip_missing_reactions:
                    continue # this reaction contains compounds that are not in the training set
                else:
                    # @@@ there is a problem for CIDs that have no InChI, I currently put nH=0, nMg=0 and z=0 for them
                    for cid in unknown_cids:
                        for (nH, z, nMg) in self.gc.cid2pseudoisomers(cid):
                            self.cid2pmap_dict.setdefault(cid, pseudoisomer.PseudoisomerMap())
                            self.cid2pmap_dict[cid].Add(nH, z, nMg, 0.0)
                        known_cids.add(cid)
            
            self.data.append(row_data)
            for cid in reaction_cids:
                self.cid2rowids.setdefault(cid, []).append(len(self.data) - 1)

        if (self.data == []):
            raise Exception("None of the measurements from NIST can be used!")
        self.train_rowids = range(len(self.data))
        self.test_rowids = range(len(self.data)) # by default assume all the training data is also the testing data

    # TODO: unfortunately, Linear Regression cannot give any initial estimation 
    # of the dG0_f because the rank of the reaction stoichiometric matrix is too low.
    # Maybe a LP approach will be more successful?
    def linear_regression(self):
        """
            Solves the dG0_f of all the compounds, by assuming there is only one dominant pseudoisomer group for
            any compound in every reaction. Selecting the right charge is done by global rules for the pKa of every
            functional group.
            The problem thus becomes a linear problem and is solved with linear regression.
        """
        
        cid_list = sorted(self.cid2pmap_dict.keys())
        cid2species = {}
        for cid in cid_list:
            comp = self.kegg.cid2compound(cid)
            cid2species[cid] = (comp.get_nH(), comp.get_charge())
        
        N = len(self.train_rowids)
        y = zeros((N, 1))
        X = zeros((N, len(cid_list)))
        for r in range(N):
            row_data = self.data[self.train_rowids[r]]
            dG0_r = row_data.dG0_r
            for (cid, coeff) in row_data.sparse.iteritems():
                c = cid_list.index(cid)
                X[r, c] = coeff
                (nH, z) = cid2species[cid]
                dG0_r -= coeff * (nH*R*row_data.T*log(10)*row_data.pH - 2.91482*(z**2 - nH)*sqrt(row_data.I) / (1 + 1.6*sqrt(row_data.I)))
            y[r, 0] = dG0_r

        inv_corr_mat = pinv(dot(X.T, X))
        dG0_f_vec = dot(dot(inv_corr_mat, X.T), y)
        
        # add the formation energies to the CID dictionary
        for c in range(len(cid_list)):
            cid = cid_list[c]
            (nH, z) = cid2species[cid]
            cid2species[cid] = (dG0_f_vec[c], nH, z)
        return cid2species

    def update_cache(self, thermodynamics, cid_to_cache=None):
        """
            update the two caches: one containing the squared error for each row in the dataset,
            and the other containing the dGf of each cid in each row (taking the stoichiometric 
            coefficient and conditions - like pH, I, T - into account)
        """
        
        rowid_list = set(self.train_rowids)
        if (cid_to_cache == None): # start over and populate the entire cache
            self.cache_cid = {}
            self.cache_error = {}
        else:
            rowid_list = rowid_list.intersection(self.cid2rowids[cid_to_cache])

        for rowid, row_data in enumerate(self.data):
            if (cid_to_cache == None):
                cid_list = row_data.sparse.keys()
            else:
                cid_list = [cid_to_cache]
                
            for cid in cid_list: 
                self.cache_cid[(rowid, cid)] = row_data.PredictFormationEnergy(thermodynamics, cid) * row_data.sparse[cid]
            dG0_pred = sum([self.cache_cid[(rowid,cid)] for cid in row_data.GetAllCids()])
            self.cache_error[rowid] = (dG0_pred - row_data.dG0_r)**2

    @staticmethod
    def sparse_reaction_to_string(sparse_reaction, kegg, cids=False):
        left = []
        right = []
        for (cid, coeff) in sparse_reaction.iteritems():
            if (cids):
                compound = "C%05d" % cid
            else:
                compound = kegg.cid2name(cid) + "(%d)" % cid
            if (coeff == 1):
                right.append(compound)
            elif (coeff > 0):
                right.append(str(coeff) + " " + compound)
            elif (coeff == -1):
                left.append(compound)
            elif (coeff < 0):
                left.append(str(-coeff) + " " + compound)
        
        return " + ".join(left) + " = " + " + ".join(right)
    
    def verify_results(self, key, thermodynamics, html_writer):
        """Calculate all the dG0_r for the reaction from NIST and compare to
           the measured data.
        
        Write results to HTML.
        
        Args:
            key: The name of this group of results.
            thermodynamics: a Thermodynamics object that provides dG estimates.
            html_writer: to write HTML.
            ignore_I: whether or not to ignore the ionic strength in NIST.
        """
        
        logging.info("calculate the correlation between %s's predictions and the NIST database" % key)
        
        known_cid_set = thermodynamics.get_all_cids()
        dG0_obs_vec = []
        dG0_est_vec = []
       
        # A mapping from each evaluation method (NIST calls separates them to
        # A, B, C and D) to the results of the relevant measurements
        evaluation_map = {}
        total_list = []
        
        cid2count = {}
        for row_data in self.data:
            for cid in row_data.GetAllCids():
                cid2count[cid] = cid2count.setdefault(cid, 0) + 1
        
        for row_data in self.data:
            unknown_set = set(row_data.GetAllCids()).difference(known_cid_set)

            if unknown_set:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.origin)
                continue
            
            #label = row_data.evaluation
            label = row_data.K_type
            
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
            n_measurements = min([cid2count[cid] for cid in row_data.GetAllCids()])
            error = abs(row_data.dG0_r - dG0_pred)

            total_list.append([error, row_data.dG0_r, dG0_pred, 
                               row_data.sparse, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               n_measurements])
        
        # plot the profile graph
        rcParams['text.usetex'] = False
        rcParams['legend.fontsize'] = 12
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = 16
        rcParams['lines.linewidth'] = 2
        rcParams['lines.markersize'] = 3
        rcParams['figure.figsize'] = [8.0, 6.0]
        rcParams['figure.dpi'] = 100
        
        fig1 = figure()
        hold(True)
        
        colors = ['purple', 'orange', 'lightgreen', 'red', 'cyan']
        for e in sorted(evaluation_map.keys()):
            (measured, predicted) = evaluation_map[e]
            label = '%s (N = %d, RMSE = %.2f [kJ/mol])' % (e, len(measured), util.calc_rmse(measured, predicted))
            c = colors.pop(0)
            plot(measured, predicted, marker='.', linestyle='None', markerfacecolor=c, markeredgecolor=c, markersize=5, label=label)
        
        legend(loc='upper left')
        
        r2 = util.calc_r2(dG0_obs_vec, dG0_est_vec)
        rmse = util.calc_rmse(dG0_obs_vec, dG0_est_vec)
        title(r'N = %d, RMSE = %.1f [kJ/mol], r$^2$ = %.2f' % (len(dG0_obs_vec), rmse, r2), fontsize=14)
        xlabel(r'$\Delta_{obs} G^\circ$ [kJ/mol]', fontsize=14)
        ylabel(r'$\Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
        min_x = min(dG0_obs_vec)
        max_x = max(dG0_obs_vec)
        plot([min_x, max_x], [min_x, max_x], 'k--')
        axis([-60, 60, -60, 60])
        
        fig2 = figure()
        hist([(row[1] - row[2]) for row in total_list], bins=arange(-50, 50, 0.5))
        rmse = util.calc_rmse(dG0_obs_vec, dG0_est_vec)
        title(r'RMSE = %.1f [kJ/mol]' % rmse, fontsize=14)
        xlabel(r'$\Delta_{obs} G^\circ - \Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
        ylabel(r'no. of measurements', fontsize=14)

        fig3 = figure()
        plot([row[9] for row in total_list], [abs(row[1] - row[2]) for row in total_list], '.')
        title(r'The effect of the number of measurements on the estimation error' % rmse, fontsize=14)
        xlabel(r'minimum no. of measurements among reaction compounds', fontsize=14)
        ylabel(r'$|| \Delta_{obs} G^\circ - \Delta_{est} G^\circ ||$ [kJ/mol]', fontsize=14)
        xscale('log')
        
        html_writer.write("<h2>%s</h2>" % key)
        
        html_writer.embed_matplotlib_figure(fig1, width=400, height=300)
        html_writer.embed_matplotlib_figure(fig2, width=400, height=300)
        
        html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (key))
        html_writer.write('<div id="%s" style="display:none">' % key)

        html_writer.embed_matplotlib_figure(fig3, width=400, height=300)

        table_headers = ["|error|", "dG0(obs)", "dG0(pred)", "reaction", "pH", "pMg", "I", "T", "evaluation", "min_num_measurements"]
        html_writer.write("<table>\n")
        html_writer.write("<tr><td>" + "</td><td>".join(table_headers) + "</td></tr>\n")
        
        for row in sorted(total_list, reverse=True):
            sparse_reaction = row[3]
            row[3] = self.kegg.sparse_to_hypertext(sparse_reaction, show_cids=False)
            html_writer.write("<tr><td>" + "</td><td>".join(["%.1f" % x for x in row[:3]] + [str(x) for x in row[3:]]) + "</td></tr>\n")
        html_writer.write("</table>\n")
        html_writer.write("</div><br>\n")

    def evaluate_MSE(self):
        return sqrt(mean(self.cache_error.values()))
        
    def leave_one_out(self, cid_to_reevaluate):
        """
            Recalculate the dG0_f (of all its species) of the compound assuming all other CIDs are known
            and using all the measurements where the selected CID appears
        """
        logging.info("reevaluation the dG0_f of C%05d ... " % cid_to_reevaluate)
        
        measurements = []
        rowid_list = set(self.train_rowids).intersection(self.cid2rowids.get(cid_to_reevaluate, []))
        if (len(rowid_list) == 0):
            raise Exception("Internal Error: there no training data involving C%05d" % cid_to_reevaluate)
        
        for rowid in rowid_list:
            (sparse_reaction, pH, I, T, unused_evaluation, dG0_r) = self.data[rowid]
            
            dG0_unknown = dG0_r # the dG of the reaction
            for cid in sparse_reaction.keys():
                if (cid != cid_to_reevaluate):
                    dG0_unknown -= self.cache_cid[(rowid,cid)] # subtract the contribution of cid to the dG0_r

            if (not cid_to_reevaluate in sparse_reaction):
                raise Exception("C%05d is not participating in the reaction on row %d" % (cid_to_reevaluate, rowid))
            dG0_unknown /= sparse_reaction[cid_to_reevaluate] # divide by the stoichiometric coefficient of the CID
            if (dG0_unknown < -10000):
                raise Exception("Internal Error: " + str((sparse_reaction, dG0_r, pH, I, T)))
            measurements.append((dG0_unknown, pH, I, T))

        old_pmap = self.cid2PseudoisomerMap(cid_to_reevaluate)
        num_pseudoisomers = len(old_pmap)
        nH_list = [nH for (nH, _z) in old_pmap.keys()]
        z_list = [z for (_nH, z) in old_pmap.keys()]
        
        try:
            dG0 = reverse_thermo.solve(measurements, array(nH_list), array(z_list))
        except reverse_thermo.ReverseTransformError as e:
            return None
        
        if (isnan(dG0).any()):
            return None
        
        new_pmap = {}
        for i in range(num_pseudoisomers):
            new_pmap[nH_list[i], z_list[i]] = [dG0[i]]
        return new_pmap
    
    def hill_climb(self, max_i):
        self.update_cache(self)
        best_MSE = self.evaluate_MSE()
        print "*** Starting point: MSE = %6.4f ***" % best_MSE
        
        # make a list of all the CIDs (except for the anchored ones) for later 
        # choosing one randomly for reevaluation
        cid_list = list( set(self.cid2rowids.keys()).difference(self.anchors) )
        
        for i in range(max_i):
            cid = cid_list[random_integers(0, len(cid_list)-1)]
            new_pmap = self.leave_one_out(cid)
            if (new_pmap != None):
                old_pmap = self.cid2PseudoisomerMap(cid)
                self.cid2pmap_dict[cid] = new_pmap
                self.update_cache(self, cid)
                curr_MSE = self.evaluate_MSE()
                if (curr_MSE > best_MSE):
                    # reject the new pmap and revert to the old one
                    self.cid2pmap_dict[cid] = old_pmap
                    self.update_cache(self, cid)
                else:
                    # update the 'best_MSE' variable
                    best_MSE = curr_MSE
            
            if (i % 1000 == 0):
                print "%05d) MSE = %6.4f" % (i, best_MSE)
        
        print "*** Ending point: MSE = %6.4f ***" % best_MSE
        
    def deterministic_hill_climb(self, max_i, verbose=True):
        """
            Unlike hill_climb, it doesn't randomly choose one CID and use all the others as a pivot to recalculate the dG0_f
            But rather tries every CID and then chooses the one with the most improvement
        """
        self.update_cache(self)
        old_MSE = self.evaluate_MSE()
        if (verbose):
            print "*** Starting point: MSE = %6.4f ***" % old_MSE
        
        # make a list of all the CIDs (except for the anchored ones) for later 
        # choosing one randomly for reevaluation
        trainable_cids = set()
        for rowid in self.train_rowids:
            sparse_reaction = self.data[rowid][0]
            for cid in sparse_reaction.keys():
                if (cid not in self.anchors):
                    trainable_cids.add(cid)
        
        for i in range(max_i):
            best_pmap = None
            best_cid = None
            best_MSE = old_MSE
            for cid in trainable_cids:
                curr_pmap = self.leave_one_out(cid)
                if (curr_pmap == None):
                    continue
                old_pmap = self.cid2PseudoisomerMap(cid)
                self.cid2pmap_dict[cid] = curr_pmap
                self.update_cache(self, cid)
                curr_MSE = self.evaluate_MSE()
                self.cid2pmap_dict[cid] = old_pmap # return to the old pmap in order to test again
                self.update_cache(self, cid)
                
                if (curr_MSE < best_MSE):
                    best_MSE = curr_MSE
                    best_cid = cid
                    best_pmap = curr_pmap
            
            if (best_pmap == None):
                logging.info("stopping because none of the CIDs could improve the MSE anymore")
                break
            self.cid2pmap_dict[best_cid] = best_pmap
            old_MSE = best_MSE
            if (verbose):
                print "%05d) C%05d - MSE = %6.4f" % (i, best_cid, old_MSE)
        
        if (verbose):
            print "*** Ending point: MSE = %6.4f ***" % old_MSE
        
    def write_pseudoisomers(self, csv_fname):
        csv_writer = csv.writer(open(csv_fname, 'w'))
        csv_writer.writerow(('cid', 'dG0_f', 'nH', 'z', 'nMg'))
        for cid in self.get_all_cids():
            for (nH, z, mgs, dG0) in self.cid2PseudoisomerMap(cid).ToMatrix:
                csv_writer.writerow((cid, dG0, nH, z, mgs))
                
    def analyse_reaction(self, reaction_to_analyze):
        known_cid_set = self.get_all_cids()
        unknown_set = set(reaction_to_analyze.keys()).difference(known_cid_set)
        if (len(unknown_set) > 0):
            raise Exception("Cannot analyze the reaction (%s) since some of the compounds is missing a dG0_f: %s " % (self.kegg.write_reaction(reaction_to_analyze), str(unknown_set)))

        # find all the rows in observed data that are the same as this reaction
        data_mat = []
        for (sparse_reaction, pH, I, T, unused_evaluation, dG0_obs) in self.test_data:
            if (reaction_to_analyze == sparse_reaction):
                dG0_est = self.reaction_to_dG0(sparse_reaction, pH, I, T)
                data_mat.append([pH, I, T, dG0_obs, dG0_est, dG0_obs-dG0_est])
        data_mat = matrix(data_mat)
        
        figure()
        rcParams['text.usetex'] = True
        rcParams['legend.fontsize'] = 6
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = 8
        rcParams['lines.linewidth'] = 0.2
        rcParams['lines.markersize'] = 2
        rc('text', usetex=True)
        rc('font', family='serif')
        plot(data_mat[:,0], data_mat[:,3], '.') # plot the dG0 vs the pH
        plot(data_mat[:,0], data_mat[:,4], '.') # plot the dG0 vs the pH
        xlabel('pH')
        ylabel('$\Delta_r G$')
        legend(['Observed', 'Estimated'])
        
        figure()
        contour(array(data_mat[:,[0,1,5]]))
        xlabel('pH')
        ylabel('I')
        
        show()

    
################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################

def main():
    db = database.SqliteDatabase('../res/gibbs.sqlite')
    html_writer = HtmlWriter("../res/nist/report.html")
    gc = GroupContribution(db)
    gc.override_gc_with_measurements = True
    gc.init()
    grad = GradientAscent(gc)
    nist = Nist(db, html_writer, gc.kegg())
    nist.FromDatabase()
    alberty = Alberty()
    hatzi = Hatzi()
    
    if True:
        grad.load_nist_data(nist, alberty, skip_missing_reactions=False, T_range=(298, 314))
        grad.verify_results("Alberty", alberty, html_writer)
        
        #grad.write_pseudoisomers("../res/nist/nist_dG0_f.csv")

        #html_writer.write("<h2>Using Group Contribution (Hatzimanikatis' implementation)</h2>")
        #html_writer.write("<h3>Correlation with the reduced NIST database (containing only compounds that appear in Alberty's list)</h3>")
        #logging.info("calculate the correlation between Hatzimanikatis' predictions and the reduced NIST database")
        #grad.verify_results("Hatzimanikatis_Reduced", hatzi, html_writer)

        #grad.load_nist_data(nist, hatzi, skip_missing_reactions=True, T_range=(298, 314))
        grad.verify_results("Hatzimanikatis", hatzi, html_writer)

        #grad.load_nist_data(nist, gc, skip_missing_reactions=True, T_range=(298, 314))
        grad.verify_results("Milo", gc, html_writer)
    elif False:
        # Run the gradient ascent algorithm, where the starting point is the same file used for training the GC algorithm
        grad.load_dG0_data("../data/thermodynamics/dG0.csv")
        # load the data for the anchors (i.e. compounds whose dG0 should not be changed - usually their value will be 0). 
        grad.anchors = grad.load_dG0_data("../data/thermodynamics/nist_anchors.csv")
        grad.load_nist_data(nist, grad, skip_missing_reactions=True)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.hill_climb(max_i=20000)
        grad.save_energies(grad.gc.comm, "gradient_cid2prm")
        grad.verify_results("gradient1")
        
    elif False:
        # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
        grad.load_nist_data(nist, alberty, skip_missing_reactions=True)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.cid2pmap_dict = alberty.cid2pmap_dict
        grad.hill_climb(max_i=20000)
        grad.save_energies(grad.gc.comm, "gradient_cid2prm")
        grad.verify_results("gradient2")
    
    elif False:
        # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
        # Use DETERMINISTIC gradient ascent
        grad.load_nist_data(nist, alberty, skip_missing_reactions=True, T_range=(24 + 273.15, 40 + 273.15))
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.cid2pmap_dict = alberty.cid2pmap_dict
        grad.deterministic_hill_climb(max_i=200)
        grad.save_energies(grad.gc.comm, "gradient_cid2prm")
        grad.verify_results("gradient_deterministic")
        
    elif False:
        # Run the gradient ascent algorithm, where the starting point arbitrary (predict all of the NIST compounds)
        grad = GradientAscent(gc)
        grad.load_nist_data(nist, skip_missing_reactions=False)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.hill_climb(max_i=20000)
        grad.save_energies(grad.gc.comm, "gradient_cid2prm")
        grad.verify_results("gradient3")
    
    elif False: # Use Alberty's table from (Mathematica 2006) to calculate the dG0 of all possible reactions in KEGG
        grad = GradientAscent(gc)
        grad.cid2pmap_dict = alberty.cid2pmap_dict
        (pH, I, T) = (7, 0, 300)
        counter = 0
        for rid in grad.kegg.get_all_rids():
            sparse_reaction = grad.kegg.rid2sparse_reaction(rid)
            try:
                dG0 = grad.reaction_to_dG0(sparse_reaction, pH, I, T, most_abundant=False)
                print "R%05d: dG0_r = %.2f [kJ/mol]" % (rid, dG0)
                counter += 1
            except MissingCompoundFormationEnergy as e:
                #print "R%05d: missing formation energy of C%05d" % (rid, e.cid)
                pass
        print "Managed to calculate the dG0 of %d reactions" % counter
        
    elif False:
        util._mkdir("../res/nist/fig")
        csv_writer = csv.writer(open("../res/nist/pseudoisomers.csv", "w"))
                
        cid_set = set()
        for row in nist.data:
            sparce_reaction = row['sparse']
            cid_set.update(sparce_reaction.keys())
        
        html_writer.write("<table border=1>\n")
        for cid in sorted(list(cid_set)):
            html_writer.write("  <tr><td>C%05d</td><td>%s</td><td>" % (cid, grad.kegg.cid2name(cid)))
            try:
                mol = grad.kegg.cid2mol(cid)
                img_fname = '../res/nist/fig/C%05d.png' % cid
                html_writer.embed_img(img_fname, "C%05d" % cid)
                mol.draw(show=False, filename=img_fname)
            except AssertionError as e:
                html_writer.write("WARNING: cannot draw C%05d - %s" % (cid, str(e)))
            except KeggParseException as e:
                html_writer.write("WARNING: cannot draw C%05d - %s" % (cid, str(e)))
            html_writer.write("</td><td>")
            if (cid in alberty.cid2pmap_dict):
                for (nH, z) in alberty.cid2pmap_dict[cid].keys():
                    html_writer.write("(nH=%d, z=%d)<br>" % (nH, z))
                    csv_writer.writerow((cid, nH, z))
            else:
                nH = grad.kegg.cid2num_hydrogens(cid)
                z = grad.kegg.cid2charge(cid)
                html_writer.write("unknown pseudoisomers<br>")
                html_writer.write("(nH=%d, z=%d)" % (nH, z))
                csv_writer.writerow((cid, nH, z))
            
            html_writer.write("</td></tr>\n")
        html_writer.write("</table>\n")
    html_writer.close()

if __name__ == '__main__':
    main()