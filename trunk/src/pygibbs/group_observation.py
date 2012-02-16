import json, logging, csv
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.group_vector import GroupVector
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.thermodynamic_constants import default_I, default_pH, default_T,\
    default_pMg, RedoxCarriers
from pygibbs.nist import Nist
from pygibbs.dissociation_constants import MissingDissociationConstantError
from toolbox.linear_regression import LinearRegression
from toolbox.html_writer import NullHtmlWriter

class GroupObservation(object):
    
    def __init__(self, obs_id, obs_type, anchored, dG0, sparse):
        self.obs_id = obs_id
        self.obs_type = obs_type # can be 'formation', 'reaction', 'acid-base' or 'Mg'
        self.dG0 = dG0
        self.sparse = sparse
        self.anchored = anchored

    def NetCharge(self):
        return self.groupvec.NetCharge()

    def Hydrogens(self):
        return self.groupvec.Hydrogens()

    def Magnesiums(self):
        return self.groupvec.Magnesiums()

    def ToDatabase(self, db, table_name):
        db.Insert(table_name, self.ToCsvRow())
    
    @staticmethod 
    def FromDatabaseRow(row):
        return GroupObservation(row['id'], row['type'], row['anchored'],
                                row['dG0'], json.loads(row['reaction']))
        
    def ToCsvRow(self):
        return [self.obs_id, self.obs_type, self.anchored, self.dG0, 
                json.dumps(self.sparse)]
    
    def __hash__(self):
        return json.dumps(self.sparse)
    
    def __eq__(self, other):
        return self.sparse == other.sparse
    
    def Normalize(self):
        """
            In self.sparse, the first non-zero coefficient must be 1,
            where the order is defined alphabetically by the keys.
        """
        factor = 1.0 / self.sparse[min(self.sparse.keys())]
        self.sparse = dict((k,v*factor) for (k,v) in self.sparse.iteritems())
        self.dG0 *= factor
        
class GroupObervationCollection(object):
    
    def __init__(self, db, html_writer, group_decomposer, dissociation, transformed=False):
        self.kegg = Kegg.getInstance()
        self.db = db
        self.groups_data = group_decomposer.groups_data
        self.group_decomposer = group_decomposer
        self.dissociation = dissociation
        self.observations = []
        self.pid2gv = {}
        self.n_groups = len(self.groups_data.all_groups)
        self.groupvec_titles = ["g%d" % i for i in xrange(self.n_groups)]
        self.transformed = transformed
        self.use_pkas = False
        self.html_writer = html_writer

        if transformed:
            prefix = 'bgc'
            self.gibbs_symbol = "&Delta;<sub>r</sub>G'&deg;"
        else:
            prefix = 'pgc'
            self.gibbs_symbol = "&Delta;<sub>r</sub>G&deg;"

        self.TRAIN_GROUPS_TABLE_NAME = prefix + '_train_groups'
        self.PSEUDOISOMER_GV_TABLE_NAME = prefix + '_pseudoisomer_groupvector'
        self.OBSERVATION_TABLE_NAME = prefix + '_observations'
    
        self.GROUP_MATRIX_TABLE_NAME = prefix + '_G'
        self.STOICHIOMETRIC_MATRIX_TABLE_NAME = prefix + '_S'
        self.GIBBS_ENERGY_TABLE_NAME = prefix + '_gibbs'
        self.PSEUDOISOMER_ID_TABLE_NAME = prefix + '_ID'
    

    
    @staticmethod
    def FromFiles(db, html_writer, group_decomposer, dissociation,
                  transformed=False, assert_decomposition=True):
        
        obs_collections = GroupObervationCollection(db, html_writer,
                                group_decomposer, dissociation, transformed)
        
        if obs_collections.use_pkas:
            html_writer.write('<br><b>List of pKa for training</b>')
            html_writer.insert_toggle(start_here=True)
            obs_collections.AddDissociationTable()
            html_writer.div_end()

        obs_collections.ReadFormationEnergies()

        html_writer.write('<br><b>List of compounds for training</b>')
        html_writer.insert_toggle(start_here=True)
        obs_collections.AddFormationEnergies()
        html_writer.div_end()

        html_writer.write('</br><b>List of NIST reactions for training</b>')
        html_writer.insert_toggle(start_here=True)
        obs_collections.AddNistDatabase(assert_decomposition=assert_decomposition)
        html_writer.div_end()
        
        obs_collections.AddRedoxCarriers()
        
        obs_collections.CalculateMatrices()
        
        return obs_collections
    
    def AddPseudoisomerFromMolecule(self, pid, mol):
        """
            Decomposes a molecule and calculates its GroupVector.
            Adds the GroupVector to the dictionary according to the given compound_id.
        """
        if pid in self.pid2gv:
            return self.pid2gv[pid]
        
        if mol is None:
            self.pid2gv[pid] = None
            return None
        
        try:
            mol.RemoveHydrogens()
            decomposition = self.group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
            groupvec = decomposition.AsVector()
            self.pid2gv[pid] = groupvec
            return groupvec
        except GroupDecompositionError:
            self.pid2gv[pid] = None
            return None

    def AddPseudoisomer(self, cid, nH, nMg):
        self.cid2nH_nMg[cid] = (nH, nMg)
        pid = self.CID2PID(cid)
        diss_table = self.GetDissociationTable(cid)
        mol = diss_table.GetMol(nH=nH, nMg=nMg)
        self.AddPseudoisomerFromMolecule(pid, mol)
        return pid
    
    def AddObservation(self, obs_id, obs_type, anchored, dG0, sparse):
        for pid in sparse.keys():
            if pid not in self.pid2gv:
                raise Exception("You must first add all relevant pseudoisomers, "
                                "and only then add the observation")
        
        obs = GroupObservation(obs_id, obs_type, anchored, dG0, sparse)
        obs.Normalize()
        self.observations.append(obs)

    def CID2PID(self, cid):
        if self.transformed:
            return "C%05d" % (cid)
        else:
            nH, nMg = self.cid2nH_nMg.get(cid, (0, 0))
            return "C%05d_nH%d_nMg%d" % (cid, nH, nMg)

    def GetDissociationTable(self, cid, assert_existing=True):
        # ToMatrix returns tuples of (nH, z, nMg, dG0)
        diss_table = self.dissociation.GetDissociationTable(cid, 
                                                create_if_missing=False)
        if assert_existing and diss_table is None:
            raise Exception("C%05d does not have a dissociation table" % (cid))
        return diss_table

    def AddDissociationTable(self):
        if self.transformed:
            raise Exception("Cannot add pKas to the training set when running"
                            " in 'transformed' mode")
        for cid in self.dissociation.GetAllCids():
            diss_table = self.GetDissociationTable(cid, assert_existing=False)
            if diss_table is None:
                continue
            name = self.kegg.cid2name(cid)
            logging.debug("Reading pKa data for C%05d, %s:" % (cid, name))
            self.html_writer.write('<b>C%05d - %s</b>\n' % (cid, name))
            for key in diss_table:
                nH_above, nH_below, nMg_above, nMg_below = key
                id_below = "C%05d_H%d_Mg%d" % (cid, nH_below, nMg_below)
                id_above = "C%05d_H%d_Mg%d" % (cid, nH_above, nMg_above)
                sparse = {id_below:1, id_above:-1}
                
                ddG0, _ref = diss_table.ddGs[key]
                mol_above = diss_table.GetMol(nH_above, nMg_above)
                mol_below = diss_table.GetMol(nH_below, nMg_below)
                
                self.AddPseudoisomerFromMolecule(id_below, mol_below)
                self.AddPseudoisomerFromMolecule(id_above, mol_above)
        
                if nH_above != nH_below:
                    s = "nH = %d -> %d" % (nH_below, nH_above)
                elif nMg_above != nMg_below:
                    s = "nMg = %d -> %d" % (nMg_below, nMg_above)
                logging.debug("\t" + s)
                html_text = ""
                html_text += '<font size="1">\n'
                html_text += '&Delta;&Delta;G = %.2f, %s</br>\n' % (ddG0, s)
                html_text += 'SMILES = %s >> %s</br>\n' % (mol_above.ToSmiles(), mol_below.ToSmiles())
                html_text += '%s >> %s</br>\n' % (mol_above.ToSVG(), mol_below.ToSVG())
                html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
                html_text += '</font>\n'
                self.html_writer.write(html_text)

                if nH_above != nH_below:
                    if nH_above != nH_below-1:
                        raise Exception("a pKa must represent a difference of exactly 1 hydrogen")
                    obs_id = '%s [nH=%d->%d]' % (name, nH_above, nH_below)
                    self.AddObservation(obs_id=obs_id, obs_type='acid-base',
                                            anchored=False, dG0=ddG0, sparse=sparse)
                elif nMg_above != nMg_below:
                    if nMg_above != nMg_below-1:
                        raise Exception("a pK_Mg must represent a difference of exactly 1 magnesium")
                    obs_id = '%s [nMg=%d->%d]' % (name, nMg_above, nMg_below)
                    raise Exception("Mg dissociation constants are currently not supported in PGC")
                    #self.AddObservation(obs_id=obs_id, obs_type='Mg',
                    #                        anchored=False, dG0=ddG0, sparse=sparse)
                    # TODO: solve the problems with Mg species !!!
                    # there seems to be a problem in decomposition of ADP:Mg (and maybe other forms that have Mg ions)

    def ReadFormationEnergies(self,
                  obs_fname='../data/thermodynamics/formation_energies.csv',
                  pH=default_pH, I=default_I, pMg=default_pMg, T=default_T):

        self.formation_dict = {}
        if self.transformed:
            self.cid2nH_nMg = {}
        else:
            self.cid2nH_nMg = self.dissociation.GetCid2nH_nMg(pH=7, I=0, pMg=14, T=default_T)
        
        for label in ['training', 'testing']:
            ptable = PsuedoisomerTableThermodynamics.FromCsvFile(obs_fname,
                                                                 label=label)
            for cid in ptable.get_all_cids():
                pmatrix = ptable.cid2PseudoisomerMap(cid).ToMatrix() 
                if len(pmatrix) != 1:
                    raise Exception("multiple training species for C%05d" % cid)
                nH, charge, nMg, dG0 = pmatrix[0]
                self.cid2nH_nMg[cid] = (nH, nMg)
                diss_table = self.GetDissociationTable(cid)
                diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
                mol = diss_table.GetMol(nH=nH, nMg=nMg)
                dG0_prime = diss_table.Transform(pH, I, pMg, T)
                ref = ptable.cid2SourceString(cid)
                self.formation_dict[cid] = (label, ref, dG0_prime,
                                            dG0, nH, charge, nMg, mol)
    
    def AddFormationEnergies(self,
            obs_fname="../data/thermodynamics/formation_energies.csv",
            pH=default_pH, I=default_I, pMg=default_pMg, T=default_T):
        """
            Add observations based on a table of derived chemical formation energies.
            If working in self.transformed mode, the dG0_f is Legendre-transformed
            to dG0'_f using the table of dissociation constants.
        """
        
        for cid, v in self.formation_dict.iteritems():
            label, ref, dG0_prime, dG0, nH, charge, nMg, mol = v
            obs_id = "%s (%d)" % (self.kegg.cid2name(cid), nH)
            logging.debug('Adding the formation energy of %s', obs_id)
            if mol is not None:
                mol.SetTitle(obs_id)
            pseudoisomer_id = self.CID2PID(cid)
            gv = self.AddPseudoisomerFromMolecule(pseudoisomer_id, mol)
            sparse = {pseudoisomer_id:1}
    
            html_text = ""
            html_text += "<b>%s [%s]</b></br>\n" % (obs_id, ref)
            html_text += '<font size="1">\n'
            html_text += 'KEGG ID = C%05d</br>\n' % cid
            if self.transformed:
                html_text += "&Delta;<sub>f</sub>G'&deg; = %.2f, " % dG0_prime
                html_text += 'pH = %g, I = %g, pMg = %g, T = %g</br>\n' % (pH, I, pMg, T)
            else:
                html_text += '&Delta;<sub>f</sub>G&deg; = %.2f, ' % dG0
                html_text += 'nH = %d, nMg = %d</br>\n' % (nH, nMg)
            if mol is not None:
                html_text += 'SMILES = %s</br>\n' % (mol.ToSmiles())
                html_text += '%s</br>\n' % (mol.ToSVG())
            html_text += 'Stoichiometry = %s</br>\n' % str(sparse)
            html_text += 'Group vector = %s</br>\n' % str(gv)

            if self.transformed:
                self.AddObservation(obs_id=obs_id, obs_type='formation',
                    anchored=(label=='testing'), dG0=dG0_prime, sparse=sparse)
            else:
                self.AddObservation(obs_id=obs_id, obs_type='formation',
                    anchored=(label=='testing'), dG0=dG0, sparse=sparse)
                
                if gv is not None:
                    # check that the nH and z of the decomposition matches
                    # the indicated values from the CSV file
                    gc_nH = gv.Hydrogens()
                    gc_charge = gv.NetCharge()
                    if nH != gc_nH:
                        s = 'ERROR: Hydrogen count doesn\'t match: explicit = %d, decomposition = %d' % (
                            nH, gc_nH)
                        logging.error(s)
                        html_text += s + '</br>\n'
                    if charge != gc_charge:
                        s = 'ERROR: Charge doesn\'t match: explicit = %d, decomposition = %d' % (
                            charge, gc_charge)
                        logging.error(s)
                        html_text += s + '</br>\n'
            html_text += '</font>\n'
            self.html_writer.write(html_text)

    def AddNistPseudoisomers(self, nist, assert_decomposition=True):
        # add compounds to the dictionary of GroupVectors
        for cid in nist.GetAllCids():
            pseudoisomer_id = self.CID2PID(cid)
            if pseudoisomer_id in self.pid2gv:
                # we already have this compound from the formation energy table
                continue

            name = self.kegg.cid2name(cid)
            nH, nMg = self.cid2nH_nMg.get(cid, (0, 0))
            if self.transformed:
                mol = self.dissociation.GetAnyMol(cid)
            else:
                mol = self.dissociation.GetMol(cid, nH, nMg)
                
            gv = self.AddPseudoisomerFromMolecule(pseudoisomer_id, mol)
            html_text = ""
            html_text += "<b>%s (%s)</b></br>\n" % (name, pseudoisomer_id)
            html_text += '<font size="1">\n'
            if mol is not None:
                html_text += 'SMILES = %s</br>\n' % (mol.ToSmiles())
                html_text += '%s</br>\n' % (mol.ToSVG())
                html_text += 'Group vector = %s</br>\n' % str(gv)
            html_text += 'Stoichiometry = %s</br>\n' % str({pseudoisomer_id:1})
            html_text += '</font>\n'
            self.html_writer.write(html_text)
        
    def AddNistDatabase(self, assert_decomposition=True):
        """
            Add the observations based on equilibrium constants from the NIST database.
            If using non-transformed group contribution, it is required to reverse Legendre-
            transform the data in order to get chemical reaction energies.
            
            This methods tries to use the same pseudoisomers for each compound as
            the table of formation energies, in order to have less columns in the 
            final stoichiometric matrix.
        """
        nist = Nist()
        self.AddNistPseudoisomers(nist, assert_decomposition)
        
        # create a dictionary from each unique reaction to the list of measured dG0'
        # and subtract from dG0' the formation energies of the anchored compounds
        for r, nist_row_data in enumerate(nist.SelectRowsFromNist()):
            obs_id = "NIST%03d" % r
            if not self.transformed:
                try:
                    dG0 = self.dissociation.ReverseTransformNistRow(
                                                nist_row_data, self.cid2nH_nMg)
                except MissingDissociationConstantError as e:
                    logging.warning('Cannot reverse transform NIST%03d because of'
                    ' of a missing dissociation constant for C%05d' % (r, e.cid))
                    continue
            else:
                dG0 = nist_row_data.dG0_r # we are using transformed energies
            
            sparse_pid = dict((self.CID2PID(cid), coeff)
                        for (cid, coeff) in nist_row_data.reaction.iteritems())

            total_gv = GroupVector(self.groups_data)
            for pseudoisomer_id, coeff in sparse_pid.iteritems():
                # ignores pseudoisomers which cannot be decomposed. These must
                # later be subtracted from the reaction in order to run PGC
                if self.pid2gv[pseudoisomer_id] is not None:
                    total_gv += self.pid2gv[pseudoisomer_id] * coeff
            
            self.AddObservation(obs_id=obs_id, obs_type='reaction',
                                    anchored=False, dG0=dG0, sparse=sparse_pid)
            
            html_text = ""
            html_text += "<b id=%s>%s</b></br>\n" % (obs_id, obs_id)
            html_text += '<font size="1">\n'
            html_text += "NIST conditions: pH = %g, I = %g, pMg = %g, T = %g</br>\n" % \
                         (nist_row_data.pH, nist_row_data.I,
                          nist_row_data.pMg, nist_row_data.T)
            html_text += 'NIST reference: <a href="%s">%s</a></br>\n' % \
                         (nist_row_data.url, nist_row_data.ref_id)
            html_text += 'EC = %s</br>\n' % nist_row_data.ec
            html_text += "Reaction: %s</br>\n" % \
                         nist_row_data.reaction.to_hypertext(show_cids=False)
            html_text += '%s: %.2f, </br>\n' % \
                         (self.gibbs_symbol, dG0)
            html_text += 'Group vector: %s</br>\n' % str(total_gv)
            html_text += '</font>\n'
            self.html_writer.write(html_text)

    def AddRedoxCarriers(self):
        redox_carriers = RedoxCarriers()
        for name, rc in redox_carriers.iteritems():
            obs_id = name + " redox"

            pid_ox = self.AddPseudoisomer(rc.cid_ox, rc.nH_ox, 0)
            pid_red = self.AddPseudoisomer(rc.cid_red, rc.nH_red, 0)
            sparse = {pid_ox:-1, pid_red:1}

            if self.transformed:
                dG0 = rc.ddG0_prime
            else:
                dG0 = rc.ddG0
            
            self.AddObservation(obs_id=obs_id, obs_type='reaction',
                                anchored=True, dG0=dG0, sparse=sparse)

    def CalculateMatrices(self, analyze_residuals=True, db_prefix=False):
        """ 
            Returns the regression matrix and corresponding observation data 
            and names as a tuple: (A, b, names)
            Next step is to solve Ax = b
            
            Makes sure that there are no duplicates (two identical rows in A)
            by averaging the observed values.
        """
        self.pid_list = sorted(self.pid2gv.keys())
        n = len(self.observations) # number of observations
        m = len(self.pid_list) # number of compounds
        g = len(self.groups_data.GetGroupNames()) # number of groups

        # first create the group matrix G (rows=pseudoisomer, cols=groups)
        self.G = np.zeros((m, g))
        for i_pid, pid in enumerate(self.pid_list):
            gv = self.pid2gv[pid]
            if gv is None:
                continue
            self.G[i_pid, :] = np.array(gv.Flatten())
                
        # then create the stoichiometric matrix S (rows=observation, cols=pseudoisomers)
        self.S = np.zeros((m, n))
        self.gibbs_values = np.zeros((1, n))
        anchored_obs = []
        unanchored_obs = []
        for i_obs, obs in enumerate(self.observations):
            for pid, coeff in obs.sparse.iteritems():
                i_pid = self.pid_list.index(pid)
                self.S[i_pid, i_obs] = coeff
            self.gibbs_values[0, i_obs] = obs.dG0
            if obs.anchored:
                anchored_obs.append(i_obs)
            else:
                unanchored_obs.append(i_obs)
                
        # now remove anchored data from S and leave only the data which will be 
        # used for calculating the group contributions
        g, _ = LinearRegression.LeastSquares(self.S[:, anchored_obs],
                                             self.gibbs_values[:, anchored_obs])
        P_C, P_L = LinearRegression.ColumnProjection(self.S[:, anchored_obs])
        anchored_gibbs_values = np.dot(np.dot(g, P_C), self.S)
        self.gibbs_values -= anchored_gibbs_values
        self.S = np.dot(P_L, self.S)
        
        # set epsilon-small values to absolute 0
        self.S[np.where(abs(self.S) < 1e-10)] = 0
        
        # removed zero rows (compounds) from S
        nonzero_rows = np.nonzero(np.sum(abs(self.S), 1))[0]
        self.S = self.S[nonzero_rows, :]
        self.G = self.G[nonzero_rows, :]
        self.pid_list = [self.pid_list[i] for i in nonzero_rows]
        
        for pid in self.pid_list:
            if self.pid2gv[pid] is None:
                raise Exception("%s has no defined structure, "
                                "but is still part of the final group matrix!")

        # removed zero column (observations) from S
        nonzero_cols = np.nonzero(np.sum(abs(self.S), 0))[0]
        self.S = self.S[:, nonzero_cols]
        self.gibbs_values = self.gibbs_values[:, nonzero_cols]
        self.observations = [self.observations[i] for i in nonzero_cols]

    def AnalyzeResiduals(self):
        GS = np.dot(self.G.T, self.S)
        # Write the analysis of residuals:
        # I am not sure if this analysis should be done before "uniquing"
        # the rows of S or after. The observation residual is much smaller
        # in the latter case, since intra-reaction noise is averaged.
        _P_R1, P_N1 = LinearRegression.RowProjection(self.S)
        _P_R2, P_N2 = LinearRegression.RowProjection(GS)
        
        r_obs = np.linalg.norm(np.dot(self.gibbs_values, P_N1))
        r_est = np.linalg.norm(np.dot(self.gibbs_values, P_N2 - P_N1))
        r_tot = np.linalg.norm(np.dot(self.gibbs_values, P_N2))
        
        self.html_writer.write('</br><b>Analysis of residuals:<b>\n')
        self.html_writer.insert_toggle(start_here=True)
        residual_text = ['r<sub>observation</sub> = %.2f kJ/mol' % r_obs,
                         'r<sub>estimation</sub> = %.2f kJ/mol' % r_est,
                         'r<sub>total</sub> = %.2f kJ/mol' % r_tot]
        self.html_writer.write_ul(residual_text)
        self.html_writer.div_end()

    def GetGroupMatrix(self):
        # 'unique' the rows S. For each set of rows that is united,
        # the Y-value for the new row is the average of the corresponding Y-values.
        GS = np.dot(self.G.T, self.S)

        unique_GS, col_mapping = LinearRegression.ColumnUnique(GS, remove_zero=True)
        unique_gibbs_values = np.zeros((1, unique_GS.shape[1]))
        unique_obs_types = []
        unique_obs_ids = []
        for i, old_indices in sorted(col_mapping.iteritems()):
            unique_gibbs_values[0, i] = np.mean(self.gibbs_values[0, old_indices])
            unique_obs_types.append(self.observations[old_indices[0]].obs_type) # take the type of the first one (not perfect...)
            unique_obs_ids.append(', '.join([self.observations[i].obs_id for i in old_indices]))            
        
        return unique_GS, unique_gibbs_values, unique_obs_ids, unique_obs_types 

    def ToDatabase(self):
        # This table is used only as an output for checking results
        # it is not used elsewhere in the code
        self.db.CreateTable(self.TRAIN_GROUPS_TABLE_NAME, 
                            'name TEXT, nH INT, z INT, nMg INT',
                            drop_if_exists=True)
        for group in self.groups_data.all_groups:
            self.db.Insert(self.TRAIN_GROUPS_TABLE_NAME,
                    [group.name, group.hydrogens, group.charge, group.nMg])

        # the PSEUDOISOMER_GV_TABLE_NAME will contain a mapping from the pseudoisomer
        # IDs to the corresponding GroupVectors
        self.db.CreateTable(self.PSEUDOISOMER_GV_TABLE_NAME,
                            'pid TEXT, groupvec TEXT',
                            drop_if_exists=True)
        for pid, gv in self.pid2gv.iteritems():
            s = None
            if gv is not None:
                s = gv.ToJSONString()
            self.db.Insert(self.PSEUDOISOMER_GV_TABLE_NAME, [pid, s])

        # the table 'group_observation' will contain all the observed data
        # that is used for training later
        self.db.CreateTable(self.OBSERVATION_TABLE_NAME,
                            'id TEXT, type TEXT, anchored BOOL, dG0 REAL, reaction TEXT',
                            drop_if_exists=True)

        for obs in self.observations:
            obs.ToDatabase(self.db, table_name=self.OBSERVATION_TABLE_NAME)
        
        self.db.SaveNumpyMatrix(self.GROUP_MATRIX_TABLE_NAME, self.G)
        # save S and gibbs_values as their transposed 
        # because otherwise they will have too many columns
        self.db.SaveNumpyMatrix(self.STOICHIOMETRIC_MATRIX_TABLE_NAME, self.S.T)
        self.db.SaveNumpyMatrix(self.GIBBS_ENERGY_TABLE_NAME, self.gibbs_values.T)

        self.db.CreateTable(self.PSEUDOISOMER_ID_TABLE_NAME, "pid TEXT")
        for pid in self.pid_list:
            self.db.Insert(self.PSEUDOISOMER_ID_TABLE_NAME, [pid])
        
        self.db.Commit()
        
    def ToCSV(self, gv_fname, obs_fname):
        # write the mapping from IDs to group vectors
        csv_writer = csv.writer(open(gv_fname, 'w'))
        csv_writer.writerow(['id', 'groupvec'])
        for pseudo_id, gv in self.pid2gv.iteritems():
            csv_writer.writerow([pseudo_id, gv.ToJSONString()])

        # write all observations
        csv_writer = csv.writer(open(obs_fname, 'w'))
        csv_writer.writerow(['id', 'type', 'anchored', 'dG0', 'reaction'])
        for obs in self.observations:
            csv_writer.writerow(obs.ToCsvRow())
        
    @staticmethod
    def FromDatabase(db, group_decomposer, transformed=False):
        html_writer = NullHtmlWriter()
        dissociation = None
        
        obs_collections = GroupObervationCollection(db, html_writer,
                                group_decomposer, dissociation, transformed)
        obs_collections.pid2gv = {}
        for row in db.DictReader(obs_collections.PSEUDOISOMER_GV_TABLE_NAME):
            if row['groupvec'] is not None:
                gv = GroupVector.FromJSONString(obs_collections.groups_data, row['groupvec'])
            else:
                gv = None
            obs_collections.pid2gv[row['pid']] = gv
        
        obs_collections.observations = []
        for row in db.DictReader(obs_collections.OBSERVATION_TABLE_NAME):
            obs = GroupObservation.FromDatabaseRow(row)
            obs_collections.observations.append(obs)

        obs_collections.G = db.LoadNumpyMatrix(obs_collections.GROUP_MATRIX_TABLE_NAME)
        obs_collections.S = db.LoadNumpyMatrix(obs_collections.STOICHIOMETRIC_MATRIX_TABLE_NAME).T
        obs_collections.gibbs_values = db.LoadNumpyMatrix(obs_collections.GIBBS_ENERGY_TABLE_NAME).T
        obs_collections.pid_list = []
        for row in db.DictReader(obs_collections.PSEUDOISOMER_ID_TABLE_NAME):
            obs_collections.pid_list.append(row['pid'])
        
        return obs_collections
        
    def ReportToHTML(self):
        n = len(self.observations) # number of observations
        m = len(self.pid_list) # number of compounds
        g = len(self.groups_data.GetGroupNames()) # number of groups

        self.html_writer.write('</br><b>Observation summary</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.html_writer.write_ul(['observations: %d' % n,
                                   'compounds: %d' % m,
                                   'groups: %d' % g])
        
        rowdicts = []
        for i in xrange(n):
            rowdict = {'#': i, 'ID': self.observations[i].obs_id}
            rowdict[self.gibbs_symbol] = self.gibbs_values[0, i]
            
            rowdict['Reaction'] = {}
            for j in np.nonzero(self.S[:, i])[0]:
                rowdict['Reaction'][self.pid_list[j]] = self.S[j, i]
            rowdicts.append(rowdict)
            
        self.html_writer.write_table(rowdicts, 
                        headers=['#', 'ID', 'Reaction', self.gibbs_symbol])
        self.html_writer.div_end()
