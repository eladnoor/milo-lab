#!/usr/bin/python

import csv
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
import pylab
import sys

from optparse import OptionParser
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.nist import Nist
from pygibbs.kegg import Kegg
from pygibbs.thermodynamics import MissingCompoundFormationEnergy, PsuedoisomerTableThermodynamics
from pygibbs.thermodynamic_constants import default_I, default_pMg, default_T,\
    default_pH, symbol_dr_G0, symbol_dr_G0_prime
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from toolbox.util import _mkdir
from pygibbs.kegg_reaction import Reaction
from matplotlib.mlab import rms_flat

class NistRegression(PsuedoisomerTableThermodynamics):
    
    def __init__(self, db, dissociation=None,
                 html_writer=None, nist=None):
        PsuedoisomerTableThermodynamics.__init__(self)
        self.db = db
        self.html_writer = html_writer or NullHtmlWriter()
        self.nist = nist or Nist()
        self.dissociation = None
        
        self.cid2pmap_dict = {}
        
        self.assume_no_pKa_by_default = False
        self.std_diff_threshold = np.inf

    def GetDissociation(self):
        if self.dissociation is None:
            self.dissociation = DissociationConstants.FromPublicDB()
        return self.dissociation

    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        
        raise MissingCompoundFormationEnergy("The compound C%05d does not have a value for its formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())
    
    def ReverseTransformWithT(self):
        """Transform preserving temperature data."""
        logging.info("Reverse transforming the NIST data")
        nist_rows = self.nist.SelectRowsFromNist()
        nist_rows_normalized = [row.Clone() for row in nist_rows]               
        data = self.GetDissociation().ReverseTranformNistRows(nist_rows_normalized)
        stoichiometric_matrix = data['S']
        cids_to_estimate = data['cids_to_estimate']
        
        logging.info("%d out of %d NIST measurements can be used" %
                     (stoichiometric_matrix.shape[0], len(nist_rows_normalized)))
        
        # for every unique row, calculate the average dG0_r of all the rows that
        # are the same reaction
        dG0_r = data['dG0_r']
        temps = data['T']
        stoich_temps = stoichiometric_matrix * temps
        stoich_and_temps = pylab.hstack((stoichiometric_matrix, stoich_temps))
        
        print stoich_and_temps.shape
        
        return stoich_and_temps, dG0_r, cids_to_estimate
        
    def ReverseTransform(self, cid2nH_nMg=None):
        """
            Performs the reverse Legendre transform on all the data in NIST where
            it is possible, i.e. where we have pKa data.
            
            Arguments:
                cid2nH_nMg - a dictionary mapping each compound ID to its chosen
                             pseudoisomer (described by nH and nMg).
        """
        logging.info("Reverse transforming the NIST data")
        nist_rows = self.nist.SelectRowsFromNist()
        logging.info("Selected %d NIST rows out of %d" %
                     (len(nist_rows), len(self.nist.data)))
        
        data = self.GetDissociation().ReverseTranformNistRows(
                                nist_rows, cid2nH_nMg=cid2nH_nMg)
                
        nist_rows_final = data['nist_rows']
        stoichiometric_matrix = data['S']
        cids_to_estimate = data['cids_to_estimate']
        n_cols = stoichiometric_matrix.shape[1]        
        
        logging.info("Only %d out of %d NIST measurements can be used" %
                     (n_cols, len(nist_rows)))

        # squeeze the regression matrix by leaving only unique rows
        unique_cols_S, col_mapping = LinearRegression.ColumnUnique(stoichiometric_matrix)
        logging.info("There are %d unique reactions" % len(col_mapping))
        unique_rids = set([nist_row.reaction.rid for nist_row in nist_rows
                            if nist_row.reaction.rid is not None])
        logging.info("Out of which %d have KEGG reaction IDs" % len(unique_rids))
        
        # for every unique column, calculate the average dG0_r of all the columns that
        # are the same reaction
        
        # full_data_mat will contain these columns: dG0, dG0_tag, dG0 - E[dG0], 
        # dG0_tag - E[dG0_tag], N
        # the averages are over the equivalence set of each reaction (i.e. the 
        # average dG of all the rows in NIST with that same reaction).
        # 'N' is the unique row number (i.e. the ID of the equivalence set)
        full_data_mat = np.matrix(np.zeros((5, n_cols)))
        full_data_mat[0, :] = np.matrix(data['dG0_r'])
        full_data_mat[1, :] = np.matrix(data['dG0_r_tag'])
        
        # unique_data_mat will contain these columns: E[dG0], E[dG0_tag],
        # std(dG0), std(dG0_tag), no. rows
        # there is exactly one row for each equivalence set (i.e. unique reaction)
        # no. rows holds the number of times this unique reaction appears in NIST
        unique_data_mat = np.matrix(np.zeros((5, len(col_mapping))))
        unique_sparse_reactions = []
        unique_nist_row_representatives = []
        for i, col_indices in col_mapping.iteritems():
            col_vector = unique_cols_S[:, i]
            
            # convert the rows of unique_rows_S to a list of sparse reactions
            sparse = {}
            for j in col_vector.nonzero()[0].flat:
                sparse[cids_to_estimate[j]] = unique_cols_S[j, i]
            reaction = Reaction(names=['NIST%03d' % i], sparse_reaction=sparse)
            unique_sparse_reactions.append(reaction)

            # find the list of indices which are equal to row i in unique_rows_S
            unique_nist_row_representatives.append(nist_rows_final[col_indices[0]])
            
            # take the mean and std of the dG0_r of these rows
            sub_data_mat  = full_data_mat[0:2, col_indices]
            unique_data_mat[0:2, i] = np.mean(sub_data_mat, 1)
            unique_data_mat[2:4, i] = np.std(sub_data_mat, 1)
            unique_data_mat[4, i]   = sub_data_mat.shape[1]
            full_data_mat[4, col_indices] = i
            full_data_mat[2:4, col_indices] = sub_data_mat
            for k in col_indices:
                # subtract the mean from each row with this reaction
                full_data_mat[2:4, k] -= unique_data_mat[0:2, i]
                    
        # write a table that lists the variances of each unique reaction
        # before and after the reverse transform
        self.WriteUniqueReactionReport(unique_sparse_reactions,
                                       unique_nist_row_representatives,
                                       unique_data_mat, full_data_mat)
        
        return unique_cols_S, unique_data_mat[0:1, :], cids_to_estimate

    @staticmethod
    def row2string(S_row, cids):
        active_cids = list(np.nonzero(S_row)[0].flat)
        sparse = dict((cids[c], S_row[c]) 
                      for c in active_cids 
                      if abs(S_row[c]) > 1e-10)
        r = Reaction("", sparse)
        return r.FullReactionString(show_cids=False)

    @staticmethod
    def row2hypertext(S_row, cids):
        kegg = Kegg.getInstance()
        active_cids = list(np.nonzero(S_row)[0].flat)
        sparse = dict((cids[c], S_row[c]) for c in active_cids)
        return kegg.sparse_to_hypertext(sparse, show_cids=False)

    def LinearRegression(self, S, obs_dG0_r, cids, cid2nH_nMg,
                         prior_thermodynamics=None):
        logging.info("Regression matrix is %d x %d" % \
                     (S.shape[0], S.shape[1]))

        cid2ref = dict((cid, 'PRC') for cid in cids)
        if prior_thermodynamics:
            # Normalize the contribution of compounds which have formation energies
            # given in the prior. Perform the regression only on the residuals
            # remaining after the normalization (note that the stoichiometric
            # matrix must also be trimmed).
            cid_index_prior = []
            dG0_prior = []
            for i, cid in enumerate(cids):
                nH, nMg = cid2nH_nMg[cid]
                try:
                    pmap_prior = prior_thermodynamics.cid2PseudoisomerMap(cid)
                except MissingCompoundFormationEnergy:
                    continue
                for p_nH, p_z, p_nMg, dG0 in pmap_prior.ToMatrix():
                    if nH == p_nH and p_nMg == nMg:
                        cid_index_prior.append(i)
                        dG0_prior.append(dG0)
                        cid2ref[cid] = pmap_prior.GetRef(p_nH, p_z, p_nMg)
                        break
            
            S_prior = np.matrix(np.zeros((len(cids), len(cid_index_prior))))
            for j, i in enumerate(cid_index_prior):
                S_prior[i, j] = 1
            dG0_prior = np.matrix(dG0_prior)
            g, _ = LinearRegression.LeastSquares(S_prior, dG0_prior)
            P_C, P_L = LinearRegression.ColumnProjection(S_prior)
            prior_dG0_r = g * P_C * S
            new_obs_dG0_r = obs_dG0_r - prior_dG0_r
            new_S = P_L * S
            
            # Find all reactions in new_S which are completely zero. This means that
            # they are completely determined by the prior.
            zero_cols = (abs(new_S).sum(0) < 1e-10).nonzero()[1]
            rowdicts = []
            for j in zero_cols.flat:
                rowdict = {}
                rowdict['reaction'] = NistRegression.row2hypertext(S[:, j], cids)
                rowdict['|error|'] = abs(new_obs_dG0_r[0, j])
                rowdict['error'] = new_obs_dG0_r[0, j]
                rowdict['NIST'] = obs_dG0_r[0, j]
                rowdict['prior'] = prior_dG0_r[0, j]
                rowdicts.append(rowdict)
            rowdicts.sort(key=lambda x:x['|error|'], reverse=True)
            self.html_writer.write('</br><b>Alberty Errors</b>\n')
            self.html_writer.write_table(rowdicts,
                                         headers=['reaction', 'error', 'NIST', 'prior'],
                                         decimal=1)
            
            est_dG0_f, _ = LinearRegression.LeastSquares(new_S, new_obs_dG0_r)
            for j, i in enumerate(cid_index_prior):
                est_dG0_f[0, i] = dG0_prior[0, j]
        else:
            est_dG0_f, _ = LinearRegression.LeastSquares(S, obs_dG0_r)
        
        est_dG0_r = est_dG0_f * S
        residuals = est_dG0_r - obs_dG0_r
        rmse = rms_flat(residuals.flat)
        logging.info("Regression results for reverse transformed data:")
        logging.info("N = %d, RMSE = %.1f" % (S.shape[1], rmse))
       
        self.html_writer.write('<p>RMSE = %.1f [kJ/mol]</p>\n' % rmse)
        rowdicts = []
        headers = ['#', 'Reaction',
                   symbol_dr_G0 + ' (obs)',
                   symbol_dr_G0 + ' (fit)',
                   symbol_dr_G0 + ' (res)']
        for i in xrange(S.shape[1]):
            rowdict = {}
            rowdict['Reaction'] = NistRegression.row2hypertext(S[:, i], cids)
            rowdict[symbol_dr_G0 + ' (obs)'] = obs_dG0_r[0, i]
            rowdict[symbol_dr_G0 + ' (fit)'] = est_dG0_r[0, i]
            rowdict[symbol_dr_G0 + ' (res)'] = residuals[0, i]
            rowdicts.append(rowdict)
        rowdicts.sort(key=lambda x:abs(x[symbol_dr_G0 + ' (res)']), reverse=True)
        self.html_writer.write_table(rowdicts, headers, decimal=1)

        # copy the solution into the diss_tables of all the compounds,
        # and then generate their PseudoisomerMaps.
        for i, cid in enumerate(cids):
            nH, nMg = cid2nH_nMg[cid]
            diss_table = self.GetDissociation().GetDissociationTable(cid)
            z = diss_table.min_charge + (nH - diss_table.min_nH)
            diss_table.SetFormationEnergyByNumHydrogens(est_dG0_f[0, i], nH, nMg)
            pmap = diss_table.GetPseudoisomerMap(nH, nMg)
            pmap.SetRef(nH, z, nMg, cid2ref[cid])
            self.cid2pmap_dict[cid] = pmap

    def WriteUniqueReactionReport(self, unique_sparse_reactions,
                                  unique_nist_row_representatives,
                                  unique_data_mat, full_data_mat,
                                  cid2nH_nMg=None):
        
        total_std = full_data_mat[2:4, :].std(1)
        
        fig = plt.figure()
        plt.plot(unique_data_mat[2, :].T, unique_data_mat[3, :].T, '.')
        plt.xlabel("$\sigma(\Delta_r G^\circ)$")
        plt.ylabel("$\sigma(\Delta_r G^{\'\circ})$")
        plt.title('$\sigma_{total}(\Delta_r G^\circ) = %.1f$ kJ/mol, '
                    '$\sigma_{total}(\Delta_r G^{\'\circ}) = %.1f$ kJ/mol' % 
                    (total_std[0, 0], total_std[1, 0]))
        self.html_writer.embed_matplotlib_figure(fig, width=640, height=480)
        logging.info('std(dG0_r) = %.1f' % total_std[0, 0])
        logging.info('std(dG\'0_r) = %.1f' % total_std[1, 0])
        
        rowdicts = []
        for i, reaction in enumerate(unique_sparse_reactions):
            logging.debug('Analyzing unique reaction: ' + 
                          str(unique_sparse_reactions[i]))
            ddG0 = self.GetDissociation().ReverseTransformReaction(reaction,
                pH=7, I=0.1, pMg=10, T=298.15, cid2nH_nMg=cid2nH_nMg)
            
            d = {}
            d["_reaction"] = reaction.to_hypertext(show_cids=False)
            d["reaction"] = reaction.FullReactionString(show_cids=False) # no hypertext for the CSV output
            d["Reference ID"] = unique_nist_row_representatives[i].ref_id
            d["EC"] = unique_nist_row_representatives[i].ec
            d["E(" + symbol_dr_G0 + ")"] = unique_data_mat[0, i]
            d["E(" + symbol_dr_G0_prime + ")"] = unique_data_mat[1, i]
            d["E(" + symbol_dr_G0 + ")'"] = unique_data_mat[0, i] + ddG0
            d["std(" + symbol_dr_G0 + ")"] = unique_data_mat[2, i]
            d["std(" + symbol_dr_G0_prime + ")"] = unique_data_mat[3, i]
            d["diff"] = unique_data_mat[2, i] - unique_data_mat[3, i]
            d["#observations"] = "%d" % unique_data_mat[4, i]
            
            flag = 0
            c_nad = reaction.sparse.get(3, 0)
            c_nadh = reaction.sparse.get(4, 0)
            c_nadp = reaction.sparse.get(6, 0)
            c_nadph = reaction.sparse.get(5, 0)
            if  c_nad == 1 and c_nadh == -1:
                flag = 1
            elif c_nad == -1 and c_nadh == 1:
                flag = -1
            elif c_nadp == 1 and c_nadph == -1:
                flag = 2
            elif c_nadp == -1 and c_nadph == 1:
                flag = -2
            d["Arren Flag"] = flag

            if d["diff"] > self.std_diff_threshold:
                _mkdir('../res/prc_reactions')
                link = "prc_reactions/%s.html" % reaction.name
                d["analysis"] = '<a href="%s">link</a>' % link
                reaction_html_writer = HtmlWriter(os.path.join('../res', link))
                self.AnalyzeSingleReaction(reaction,
                                           html_writer=reaction_html_writer)
            rowdicts.append(d)
        
        result_headers = ["E(" + symbol_dr_G0 + ")",
                          "E(" + symbol_dr_G0_prime + ")", 
                          "E(" + symbol_dr_G0 + ")'",
                          "std(" + symbol_dr_G0 + ")",
                          "std(" + symbol_dr_G0_prime + ")"]
        rowdicts.sort(key=lambda x:x["diff"], reverse=True)
        self.html_writer.write_table(rowdicts, ["reaction", "Reference ID"] + 
                                     result_headers + ["EC", "#observations", "analysis"],
                                     decimal=1)
        csv_writer = csv.DictWriter(open('../res/nist_regression_unique.csv', 'w'),
                                    ["_reaction", "Reference ID", "EC", "#observations"]
                                    + result_headers + ['Arren Flag'],
                                    extrasaction='ignore')
        csv_writer.writeheader()
        csv_writer.writerows(rowdicts)

    def AnalyzeSingleReaction(self, reaction, html_writer=None):
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.size'] = 8
        plt.rcParams['lines.linewidth'] = 2
        plt.rcParams['lines.markersize'] = 8

        if not html_writer:
            html_writer = self.html_writer

        # gather all the measurements from NIST that correspond to this reaction
        logging.info("Analyzing reaction - " + str(reaction))
        nist_rows = self.nist.SelectRowsFromNist(reaction)
        logging.info("Found %d measurements in NIST with this reaction" %
                     len(nist_rows))
        html_writer.write('<h2>%s</h2>\n' % reaction.name)
        html_writer.write('<p>\nShow observation table: ')
        div_id = html_writer.insert_toggle()
        html_writer.div_start(div_id)
        dict_list = []
        for nist_row_data in nist_rows:
            d = {}
            d['pH'] = nist_row_data.pH
            d['I'] = nist_row_data.I
            d['pMg'] = nist_row_data.pMg
            d['dG\'0_r'] = "%.2f" % nist_row_data.dG0_r
            d['T(K)'] = nist_row_data.T
            if nist_row_data.url:
                d['URL'] = '<a href="%s">link</a>' % nist_row_data.url
            else:
                d['URL'] = ''
            dict_list.append(d)
        html_writer.write_table(dict_list, headers=['T(K)', 'pH', 'I', 'pMg', 'dG\'0_r', 'URL'])
        html_writer.div_end()
        html_writer.write('</p>\n')

        # reverse transform the data
        data = self.GetDissociation().ReverseTranformNistRows(nist_rows)
        
        html_writer.write('Reaction: %s</br>\n' % \
                          reaction.to_hypertext(show_cids=False))
        fig1 = plt.figure(figsize=[9.0, 6.0], dpi=80)
        html_writer.write('Standard deviations:</br>\n<ul>\n')
        
        y_labels = ['$\Delta_r G^{\'\circ}$', '$\Delta_r G^{\circ}$']
        x_limits = {'pH' : (3, 12), 'I' : (0, 1), 'pMg' : (0, 14)}
        
        for j, y_axis in enumerate(['dG0_r_tag', 'dG0_r']):
            sigma = np.std(data[y_axis])
            html_writer.write("  <li>stdev(%s) = %.2g</li>" % (y_axis, sigma))
            for i, x_axis in enumerate(['pH', 'I', 'pMg']):
                plt.subplot(2,3,i+3*j+1)
                plt.plot(data[x_axis], data[y_axis], '.')
                plt.xlim(x_limits[x_axis])
                if j == 1:
                    plt.xlabel(x_axis)
                if i == 0:
                    plt.ylabel(y_labels[j])
        html_writer.write('</ul>\n')
        html_writer.embed_matplotlib_figure(fig1)
        
        # draw the response of the graph to pH, I and pMg:
        fig2 = plt.figure(figsize=[9.0, 5.0], dpi=80)

        pH_range = np.arange(3, 12.01, 0.05)
        I_range = np.arange(0.0, 1.01, 0.01)
        pMg_range = np.arange(0.0, 14.01, 0.1)

        ddG_vs_pH = []
        for pH in pH_range:
            ddG = self.GetDissociation().ReverseTransformReaction(reaction, pH=pH, I=default_I, 
                                                pMg=default_pMg, T=default_T)
            ddG_vs_pH.append(ddG)
        
        ddG_vs_I = []
        for I in I_range:
            ddG = self.GetDissociation().ReverseTransformReaction(reaction, pH=default_pH, I=I, 
                                                pMg=default_pMg, T=default_T)
            ddG_vs_I.append(ddG)

        ddG_vs_pMg = []
        for pMg in pMg_range:
            ddG = self.GetDissociation().ReverseTransformReaction(reaction, pH=default_pH, I=default_I, 
                                                pMg=pMg, T=default_T)
            ddG_vs_pMg.append(ddG)
        
        plt.subplot(1, 3, 1)
        plt.plot(pH_range, ddG_vs_pH, 'g-')
        plt.xlabel('pH')
        plt.ylabel('$\Delta_r G^{\'\circ} - \Delta_r G^{\circ}$')
        plt.subplot(1, 3, 2)
        plt.plot(I_range, ddG_vs_I, 'b-')
        plt.xlabel('I')
        plt.subplot(1, 3, 3)
        plt.plot(pMg_range, ddG_vs_pMg, 'r-')
        plt.xlabel('pMg')
        html_writer.embed_matplotlib_figure(fig2, width=640, height=480)
        html_writer.write('</br>\n')
        
    def ConvertPseudoisomer(self, cid, dG0, nH_from, nH_to=None):
        return self.GetDissociation().ConvertPseudoisomer(cid, dG0, nH_from, nH_to)
    
    def VerifyResults(self):
        return self.nist.verify_results(html_writer=self.html_writer, 
                                        thermodynamics=self)

    def Train(self, FromDatabase=True, prior_thermodynamics=None):
        if FromDatabase and self.db.DoesTableExist('prc_S'):
            S = self.db.LoadSparseNumpyMatrix('prc_S')
            dG0 = self.db.LoadNumpyMatrix('prc_b').T
            cids = []
            cid2nH_nMg = {}
            for rowdict in self.db.DictReader('prc_compounds'):
                cid, nH, nMg = int(rowdict['cid']), int(rowdict['nH']), int(rowdict['nMg'])
                cids.append(int(rowdict['cid']))
                cid2nH_nMg[cid] = (nH, nMg)
        else:
            cid2nH_nMg = self.GetDissociation().GetCid2nH_nMg(
                                            self.pH, self.I, self.pMg, self.T)
            S, dG0, cids = self.ReverseTransform(cid2nH_nMg=cid2nH_nMg)
            self.db.SaveSparseNumpyMatrix('prc_S', S)
            self.db.SaveNumpyMatrix('prc_b', dG0.T)
            self.db.CreateTable('prc_compounds',
                                'cid INT, name TEXT, nH INT, nMg INT')
            kegg = Kegg.getInstance()
            for cid in cids:
                nH, nMg = cid2nH_nMg[cid]
                self.db.Insert('prc_compounds',
                               [cid, kegg.cid2name(cid), nH, nMg])
            self.db.Commit()

        # Train the formation energies using linear regression
        self.LinearRegression(S, dG0, cids, cid2nH_nMg, prior_thermodynamics)
        self.ToDatabase(self.db, 'prc_pseudoisomers')
    
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-d", "--from_database", action="store_true",
                          dest="from_database", default=False,
                          help="A flag for loading the data from the cache DB instead of "
                               "the NIST directly (saves the time it takes to reverse transform).")
    opt_parser.add_option("-p", "--use_prior", action='store_true',
                          dest="use_prior", default=False,
                          help="The name of the thermodynamics file to load.")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          default='../res/prc.html',
                          help="Where to write output to.")
    return opt_parser


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    
    db = SqliteDatabase("../res/gibbs.sqlite")
    public_db = SqliteDatabase("../data/public_data.sqlite")
    output_filename = os.path.abspath(options.output_filename)
    logging.info('Will write output to %s' % output_filename)
    
    html_writer = HtmlWriter(output_filename)
    nist = Nist(T_range=None)
    nist_regression = NistRegression(db, html_writer=html_writer, nist=nist)
    nist_regression.std_diff_threshold = 5 # the threshold over which to print an analysis of a reaction
    #nist_regression.nist.T_range = None(273.15 + 24, 273.15 + 40)
    #nist_regression.nist.override_I = 0.25
    #nist_regression.nist.override_pMg = 14.0

    html_writer.write("<h2>NIST regression:</h2>")
    if options.use_prior:
        logging.info('Using the data from Alberty as fixed prior')
        prior_thermo = PsuedoisomerTableThermodynamics.FromDatabase(
            public_db, 'alberty_pseudoisomers', name="Alberty")
    else:
        prior_thermo = None
    html_writer.write('</br><b>Regression Tables</b>\n')
    html_writer.insert_toggle(start_here=True)
    nist_regression.Train(options.from_database, prior_thermo)
    html_writer.div_end()
 
    html_writer.write('</br><b>PRC results</b>\n')
    html_writer.insert_toggle(start_here=True)
    nist_regression.WriteDataToHtml(html_writer)
    html_writer.div_end()

    html_writer.write('</br><b>Transformed reaction energies - PRC vs. Observed</b>\n')
    html_writer.insert_toggle(start_here=True)
    N, rmse = nist_regression.VerifyResults()
    html_writer.div_end()
    
    logging.info("Regression results for transformed data:")
    logging.info("N = %d, RMSE = %.1f" % (N, rmse))

    html_writer.close()
    
if __name__ == "__main__":
    main()
