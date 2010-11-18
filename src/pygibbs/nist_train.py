import csv
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from toolbox.html_writer import HtmlWriter
from pygibbs.groups import GroupContribution
from pygibbs.kegg import KeggParseException
from pygibbs.thermodynamics import R,  default_I, Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.alberty import Alberty
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.nist import Nist
from pygibbs import reverse_thermo
from toolbox import util

class GradientAscent (Thermodynamics):
    def __init__(self, gc):
        Thermodynamics.__init__(self)
        self.source_string = "Gradient Ascent"
        self.gc = gc
        self.kegg = gc.kegg()
        self.data = []           # this will hold all the training data in one place, each row represents a measurement
        self.train_rowids = []  # the row indices in self.data to use for training
        self.test_rowids = []   # the row indices in self.data to use for testing
        self.cid2rowids = {} # holds a list of indices to self.train_data where the CID is participating
        self.cid2pmap_dict = {}       # this is the solution vector (what we are trying to solve)
        self.anchors = set()     # the CIDs of the anchor compounds (i.e. their pmap must be fixed)
        self.cache_error = {}    # for each rowid - holds the last calculated squared error (the difference between calc and measured)
        self.cache_cid = {}      # for each (rowid,cid) pair - holds the last calculated dG_f (multiplied by the stoichiometric coeff)

    def cid2pmap(self, cid):
        if (cid in self.cid2pmap_dict):
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
        csv_reader = csv.reader(open(train_csv_fname))
        titles = csv_reader.next() # "cid","compound name","dG0","charge","hydrogens"
        try:
            i_cid = titles.index("cid")
            i_dG0 = titles.index("dG0")
            i_z = titles.index("charge")
            i_nH = titles.index("hydrogens")
        except ValueError:
            raise Exception("this CSV file %s does not have the correct titles" % train_csv_fname)
        
        try:
            i_use = titles.index("use_for")
        except ValueError:
            i_use = None
        
        cids = set()
        for row in csv_reader:
            if (format == 'long'):
                if (i_use != None and row[i_use] in ['skip']):
                    continue
                cid = int(row[i_cid])
                dG0 = float(row[i_dG0])
                z = int(row[i_z])
                nH = int(row[i_nH])
            
            if (cid > 0):
                self.cid2pmap_dict.setdefault(cid, {})[nH, z] = [dG0]
                cids.add(cid)
        return cids

    def load_nist_data(self, nist, thermodynamics, override_I=None, skip_missing_reactions=False, T_range=None):
        """
            skip_missing_reactions - if True, the method will not add reactions which have a compound
            with a no data (on dG0_f). In this case, the GradientAscent object is meant only to improve
            the dG0_f of the training data, but not to discover formation energies of new compounds.
        """
        self.data = []
        self.cid2rowids = {}
        
        known_cids = set(thermodynamics.get_all_cids())
        for row in nist.data:
            sparse_reaction = row[6]
            [Keq, T, I, pH] = row[8:12]
            if (T_range != None and not (T_range[0] < T < T_range[1])):
                continue # the temperature is outside the allowed range
            dG0_r = -R*T*log(Keq)
            if (override_I != None):
                I = override_I
            if (row[3] == 'A'):
                evaluation = 'A'
            else:
                evaluation = 'B - D'
            #evaluation = row[2]
            reaction_cids = set(sparse_reaction.keys())
            unknown_cids = reaction_cids.difference(known_cids)
            
            if len(unknown_cids) > 0:
                if skip_missing_reactions:
                    continue # this reaction contains compounds that are not in the training set
                else:
                    for cid in unknown_cids: # @@@ there is a problem for CIDs that have no InChI, I currently put nH=0, z=0 for them
                        for (nH, z) in self.gc.cid2pseudoisomers(cid):
                            self.cid2pmap_dict.setdefault(cid, {})[nH, z] = [0.0]
                        known_cids.add(cid)
            
            self.data.append((sparse_reaction, pH, I, T, evaluation, dG0_r))
            rowid = len(self.data) - 1
            for cid in reaction_cids:
                self.cid2rowids.setdefault(cid, []).append(rowid)

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
            (sparse_reaction, pH, I, T, evaluation, dG0_r_transformed) = self.data[self.train_rowids[r]]
            dG0_r = dG0_r_transformed
            for (cid, coeff) in sparse_reaction.iteritems():
                c = cid_list.index(cid)
                X[r, c] = coeff
                (nH, z) = cid2species[cid]
                dG0_r -= coeff * (nH*R*T*log(10)*pH - 2.91482*(z**2 - nH)*sqrt(I) / (1 + 1.6*sqrt(I)))
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

        for rowid in rowid_list:
            (sparse_reaction, pH, I, T, evaluation, dG0_r) = self.data[rowid]
            if (cid_to_cache == None):
                cid_list = sparse_reaction.keys()
            else:
                cid_list = [cid_to_cache]
                
            for cid in cid_list: 
                self.cache_cid[(rowid,cid)] = Thermodynamics.pmap_to_dG0(thermodynamics.cid2pmap(cid), pH, I, T) * sparse_reaction[cid]
            dG0_pred = sum([self.cache_cid[(rowid,cid)] for cid in sparse_reaction.keys()])
            self.cache_error[rowid] = (dG0_pred - dG0_r)**2

    def save_energies(self):
        self.gc.comm.execute("DROP TABLE IF EXISTS cid2prm")
        self.gc.comm.execute("CREATE TABLE cid2prm (cid INT, dG0_f REAL, nH INT, z INT, anchor BOOL)")
        for cid in self.get_all_cids():
            for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(self.cid2pmap(cid)):
                self.gc.comm.execute("INSERT INTO cid2prm VALUES(?,?,?,?,?)", (cid, dG0, nH, z, cid in self.anchors))
        self.gc.comm.commit()

    def load_energies(self):
        self.cid2pmap_dict = {}
        self.anchors = set()
        for row in self.gc.comm.execute("SELECT * FROM cid2prm"):
            (cid, dG0, nH, z, anchor) = row
            self.cid2pmap_dict.setdefault(cid, {})
            self.cid2pmap_dict[cid].setdefault(nH, z, []).append(dG0)
            if (anchor):
                self.anchors.add(cid)
        self.update_cache(self)

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
    
    def verify_results(self, key, thermodynamics, html_writer, ignore_I=False):
        """
            recalculate all the dG0_r for the reaction from NIST and compare to the measured data
        """
        
        sys.stderr.write("Calculate the correlation between %s's predictions and the NIST database\n" % key)
        
        known_cid_set = thermodynamics.get_all_cids()
        dG0_obs_vec = []
        dG0_est_vec = []
        evaluation_map = {} # a mapping from each evaluation method (NIST calls separates them to A, B, C and D) to the results of the relevant measurements
        total_list = []
        
        cid2count = {}
        for rowid in self.train_rowids:
            (sparse_reaction, pH, I, T, evaluation, dG0_r) = self.data[rowid]
            for cid in sparse_reaction.keys():
                cid2count[cid] = cid2count.setdefault(cid, 0) + 1
        
        for rowid in self.test_rowids:
            (sparse_reaction, pH, I, T, evaluation, dG0_r) = self.data[rowid]
            if (ignore_I):
                I = default_I
            unknown_set = set(sparse_reaction.keys()).difference(known_cid_set)
            if (len(unknown_set) > 0):
                sys.stderr.write("One of the compounds in reaction at row %d in NIST doesn't have a dG0_f\n" % rowid)
            if (evaluation not in evaluation_map):
                evaluation_map[evaluation] = ([], [])
            
            try:
                dG0_pred = thermodynamics.reaction_to_dG0(sparse_reaction, pH, I, T)
            except MissingCompoundFormationEnergy:
                sys.stderr.write("One of the compounds in reaction at row %d in NIST doesn't have a dG0_f\n" % rowid)
                continue
                
            dG0_obs_vec.append(dG0_r)
            dG0_est_vec.append(dG0_pred)
            evaluation_map[evaluation][0].append(dG0_r)
            evaluation_map[evaluation][1].append(dG0_pred)
            n_measurements = min([cid2count[cid] for cid in sparse_reaction.keys()])
            total_list.append([abs(dG0_r - dG0_pred), dG0_r, dG0_pred, sparse_reaction, pH, I, T, evaluation, n_measurements])
        
        fig1 = figure()
        hold(True)
        leg = []
        
        for e in sorted(evaluation_map.keys()):
            (measured, predicted) = evaluation_map[e]
            plot(measured, predicted, '+')
            leg.append('%s (N = %d, RMSE = %.2f [kJ/mol])' % (e, len(measured), util.calc_rmse(measured, predicted)))
        
        prop = font_manager.FontProperties(size=12)
        legend(leg, loc='upper left', prop=prop)
        
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
        plot([row[8] for row in total_list], [abs(row[1] - row[2]) for row in total_list], '.')
        title(r'The effect of the number of measurements on the estimation error' % rmse, fontsize=14)
        xlabel(r'minimum no. of measurements among reaction compounds', fontsize=14)
        ylabel(r'$|| \Delta_{obs} G^\circ - \Delta_{est} G^\circ ||$ [kJ/mol]', fontsize=14)
        xscale('log')
        
        html_writer.write("<h2>%s</h2>" % key)
        
        html_writer.embed_matplotlib_figure(fig1, width=640, height=480)
        html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
        
        html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (key))
        html_writer.write('<div id="%s" style="display:none">' % key)

        html_writer.embed_matplotlib_figure(fig3, width=640, height=480)

        table_headers = ["|error|", "dG0(obs)", "dG0(pred)", "reaction", "pH", "I", "T", "evaluation", "min_num_measurements"]
        html_writer.write("<table>\n")
        html_writer.write("<tr><td>" + "</td><td>".join(table_headers) + "</td></tr>\n")
        
        for row in sorted(total_list, reverse=True):
            sparse_reaction = row[3]
            row[3] = self.kegg.sparse_to_hypertext(sparse_reaction, show_cids=False)
            html_writer.write("<tr><td>" + "</td><td>".join([str(x) for x in row]) + "</td></tr>\n")
        html_writer.write("</table>\n")
        html_writer.write("</div><br>\n")

    def evaluate_MSE(self):
        return sqrt(mean(self.cache_error.values()))
        
    def leave_one_out(self, cid_to_reevaluate):
        """
            Recalculate the dG0_f (of all its species) of the compound assuming all other CIDs are known
            and using all the measurements where the selected CID appears
        """
        #sys.stderr.write("Reevaluation the dG0_f of C%05d ... " % cid_to_reevaluate)
        
        measurements = []
        rowid_list = set(self.train_rowids).intersection(self.cid2rowids.get(cid_to_reevaluate, []))
        if (len(rowid_list) == 0):
            raise Exception("Internal Error: there no training data involving C%05d" % cid_to_reevaluate)
        
        for rowid in rowid_list:
            (sparse_reaction, pH, I, T, evaluation, dG0_r) = self.data[rowid]
            
            dG0_unknown = dG0_r # the dG of the reaction
            for (cid, coeff) in sparse_reaction.iteritems():
                if (cid != cid_to_reevaluate):
                    dG0_unknown -= self.cache_cid[(rowid,cid)] # subtract the contribution of cid to the dG0_r

            if (not cid_to_reevaluate in sparse_reaction):
                raise Exception("C%05d is not participating in the reaction on row %d" % (cid_to_reevaluate, rowid))
            dG0_unknown /= sparse_reaction[cid_to_reevaluate] # divide by the stoichiometric coefficient of the CID
            if (dG0_unknown < -10000):
                raise Exception("Internal Error: " + str((sparse_reaction, dG0_r, pH, I, T)))
            measurements.append((dG0_unknown, pH, I, T))

        old_pmap = self.cid2pmap(cid_to_reevaluate)
        num_pseudoisomers = len(old_pmap)
        nH_list = [nH for (nH, z) in old_pmap.keys()]
        z_list = [z for (nH, z) in old_pmap.keys()]
        
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
                old_pmap = self.cid2pmap(cid)
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
                old_pmap = self.cid2pmap(cid)
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
                sys.stderr.write(" Stopping because none of the CIDs could improve the MSE anymore")
                break
            self.cid2pmap_dict[best_cid] = best_pmap
            old_MSE = best_MSE
            if (verbose):
                print "%05d) C%05d - MSE = %6.4f" % (i, best_cid, old_MSE)
        
        if (verbose):
            print "*** Ending point: MSE = %6.4f ***" % old_MSE
        
    def write_pseudoisomers(self, csv_fname):
        csv_writer = csv.writer(open(csv_fname, 'w'))
        csv_writer.writerow(('CID', 'dG0_f', 'nH', 'charge'))
        for cid in self.get_all_cids():
            for (nH, z, dG0) in Thermodynamics.pmap_to_matrix(self.cid2pmap(cid)):
                csv_writer.writerow((cid, dG0, nH, z))
                
    def analyse_reaction(self, reaction_to_analyze):
        known_cid_set = self.get_all_cids()
        unknown_set = set(reaction_to_analyze.keys()).difference(known_cid_set)
        if (len(unknown_set) > 0):
            raise Exception("Cannot analyze the reaction (%s) since some of the compounds is missing a dG0_f: %s " % (self.kegg.write_reaction(reaction_to_analyze), str(unknown_set)))

        # find all the rows in observed data that are the same as this reaction
        data_mat = []
        for (sparse_reaction, pH, I, T, evaluation, dG0_obs) in self.test_data:
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

    gc = GroupContribution(sqlite_name="gibbs.sqlite")
    gc.init()
    grad = GradientAscent(gc)
    nist = Nist(gc.kegg())
    alberty = Alberty()
    hatzi = Hatzi()
    
    if True:
        html_writer = HtmlWriter("../res/nist/report.html")

        grad.load_nist_data(nist, alberty, skip_missing_reactions=False, T_range=(298, 314))
        grad.verify_results("Alberty", alberty, html_writer)
        
        #grad.write_pseudoisomers("../res/nist/nist_dG0_f.csv")

        #html_writer.write("<h2>Using Group Contribution (Hatzimanikatis' implementation)</h2>")
        #html_writer.write("<h3>Correlation with the reduced NIST database (containing only compounds that appear in Alberty's list)</h3>")
        #sys.stderr.write("Calculate the correlation between Hatzimanikatis' predictions and the reduced NIST database\n")
        #grad.verify_results("Hatzimanikatis_Reduced", hatzi, html_writer)

        #grad.load_nist_data(nist, hatzi, skip_missing_reactions=True, T_range=(298, 314))
        grad.verify_results("Hatzimanikatis", hatzi, html_writer)

        #grad.load_nist_data(nist, gc, skip_missing_reactions=True, T_range=(298, 314))
        grad.verify_results("Milo", gc, html_writer)
        html_writer.close()
        
    elif False:
        # Run the gradient ascent algorithm, where the starting point is the same file used for training the GC algorithm
        grad.load_dG0_data("../data/thermodynamics/dG0.csv")
        # load the data for the anchors (i.e. compounds whose dG0 should not be changed - usually their value will be 0). 
        grad.anchors = grad.load_dG0_data("../data/thermodynamics/nist_anchors.csv")
        grad.load_nist_data(nist, grad, skip_missing_reactions=True)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.hill_climb(max_i=20000)
        grad.save_energies()
        grad.verify_results("gradient1")
        
    elif False:
        # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
        grad.load_nist_data(nist, alberty, skip_missing_reactions=True)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.cid2pmap_dict = alberty.cid2pmap_dict
        grad.hill_climb(max_i=20000)
        grad.save_energies()
        grad.verify_results("gradient2")
    
    elif False:
        # Run the gradient ascent algorithm, where the starting point is Alberty's table from (Mathematica 2006)
        # Use DETERMINISTIC gradient ascent
        grad.load_nist_data(nist, alberty, skip_missing_reactions=True, T_range=(24 + 273.15, 40 + 273.15))
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.cid2pmap_dict = alberty.cid2pmap_dict
        grad.deterministic_hill_climb(max_i=200)
        grad.save_energies()
        grad.verify_results("gradient_deterministic")
        
    elif False:
        # Run the gradient ascent algorithm, where the starting point arbitrary (predict all of the NIST compounds)
        grad = GradientAscent(gc)
        grad.load_nist_data(nist, skip_missing_reactions=False)
        print "Training %d compounds using %d reactions: " % (len(grad.cid2pmap_dict.keys()), len(grad.data))
        grad.hill_climb(max_i=20000)
        grad.save_energies()
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
        html_writer = HtmlWriter("../res/nist/pseudoisomers.html")
        csv_writer = csv.writer(open("../res/nist/pseudoisomers.csv", "w"))
                
        cid_set = set()
        for row in nist.data[1:]:
            sparce_reaction = row[6]
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

if (__name__ == '__main__'):
    main()