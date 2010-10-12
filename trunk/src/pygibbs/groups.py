import sys, pybel, openbabel, csv, os, pylab, sqlite3, re
from copy import deepcopy
from matplotlib.font_manager import FontProperties
from toolbox.html_writer import HtmlWriter
from toolbox.util import matrixrank, multi_distribute, log_sum_exp
from pygibbs.thermodynamics import R, default_pH, default_I, default_T, default_c0, Thermodynamics
from pygibbs.feasibility import find_feasible_concentrations, LinProgNoSolutionException, find_pCr
from pygibbs import kegg

def find_smarts(smarts_str, mol):
    """
        This corrects the pyBel version of Smarts.findall() which returns results as tuples,
        and as 1-based indices, although Molecule.atoms is 0-based.
        Note that 'nodes' is 1-based, since that's how Smarts works
    """
    results = []
    for match in pybel.Smarts(smarts_str).findall(mol):
        results.append([(n - 1) for n in match])
    return results

def find_pchains(mol, lengths=[1, 2, 3], ignore_protonations=False):
    """
        end should be 'OC' for chains that do not really end, but link to carbons
        end should be '[O-1,OH]' for chains that end in an hydroxyl
    """
    group_map = {}
    group_map[("-OPO3-", 1, 0)] = [] # init group
    group_map[("-OPO3-", 0, -1)] = [] # init group
    group_map[("-OPO2-", 1, 0)] = [] # middle group
    group_map[("-OPO2-", 0, -1)] = [] # middle group
    #group_map[("-OPO3",  2,  0)] = [] # final group
    group_map[("-OPO3", 1, -1)] = [] # final group
    group_map[("-OPO3", 0, -2)] = [] # final group

    v_charge = [a.formalcharge for a in mol.atoms]
    
    for length in lengths:
        smarts_str = "CO" + ("P(=O)([OH,O-])O" * length) + "C"
        for pchain in find_smarts(smarts_str, mol):
            if (ignore_protonations):
                group_map[("-OPO3-", 0, -1)].append(set(pchain[1:6]))
            else:
                charge = v_charge[pchain[4]]
                protons = charge + 1
                group_map[("-OPO3-", protons, charge)].append(set(pchain[1:6]))
            for i in range(6, len(pchain) - 1, 4):
                if (ignore_protonations):
                    group_map[("-OPO2-", 0, -1)].append(set(pchain[i:(i + 4)]))
                else:
                    charge = v_charge[pchain[i + 2]]
                    protons = charge + 1
                    group_map[("-OPO2-", protons, charge)].append(set(pchain[i:(i + 4)]))
    
        smarts_str = "[OH,O-]" + ("P(=O)([OH,O-])O" * length) + "C"
        for pchain in find_smarts(smarts_str, mol):
            if (ignore_protonations):
                group_map[("-OPO3", 1, -1)].append(set(pchain[0:5]))
            else:
                charge = v_charge[pchain[0]] + v_charge[pchain[3]]
                protons = charge + 2
                if (("-OPO3", protons, charge) not in group_map):
                    sys.stderr.write("WARNING: This protonation (%d) level is not allowed for terminal phosphate groups. " % protons)
                    sys.stderr.write("Assuming the level is actually 1 - i.e. the charge is (-1).\n")
                    group_map[("-OPO3", 1, -1)].append(set(pchain[0:5]))
                else:
                    group_map[("-OPO3", protons, charge)].append(set(pchain[0:5]))
            for i in range(5, len(pchain) - 1, 4):
                if (ignore_protonations):
                    group_map[("-OPO2-", 0, -1)].append(set(pchain[i:(i + 4)]))
                else:
                    charge = v_charge[pchain[i + 2]]
                    protons = charge + 1
                    group_map[("-OPO2-", protons, charge)].append(set(pchain[i:(i + 4)]))
                
    return sorted(list(group_map.iteritems()))

class GroupDecompositionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class GroupMissingTrainDataError(Exception):
    def __init__(self, value, missing_groups=[]):
        self.value = value
        self.missing_groups = missing_groups
    def __str__(self):
        return repr(self.value)

class GroupContribution:    
    def __init__(self, sqlite_name, html_name, log_file=sys.stderr):
        self._kegg = None
        self.LOG_FILE = log_file
        if (not os.path.exists("../res")):
            os.mkdir("../res")
        
        self.EXP_NAME = html_name
        self.HTML = HtmlWriter("../res/" + html_name + ".html")
        self.FIG_DIR = '../res/' + self.EXP_NAME
        if (not os.path.exists(self.FIG_DIR)):
            os.mkdir(self.FIG_DIR)
        
        self.comm = sqlite3.connect("../res/" + sqlite_name)
        self.use_measured_train_values = True
    
    def __del__(self):
        self.comm.close()
        self.HTML.close()
    
    def write_gc_tables(self):
        table_names = ["groups", "contribution", "observation"]
        self.HTML.write('<ul>\n')
        for table_name in table_names:
            self.HTML.write('  <li><a href="#%s">Table %s from the database</a></li>\n' % (table_name, table_name))
        self.HTML.write('</ul>\n')
        for table_name in table_names:
            self.write_table_to_html(self.HTML, table_name)        
    
    def init(self):
        self.load_groups()
        self.load_contributions()
        self.load_concentrations()
        self.load_training_data()
        self.load_hatzimanikatis_rid_data()
        self.load_hatzimanikatis_cid_data()
        self.load_cid2pmap()

    def load_cid2pmap(self):
        self.cid2pmap = {}
        if (self.does_table_exist('gc_cid2prm')):
            for row in self.comm.execute("SELECT cid, nH, z, dG0 from gc_cid2prm;"):
                (cid, nH, z, dG0) = row
                if (cid not in self.cid2pmap):
                    self.cid2pmap[cid] = {}
                self.cid2pmap[cid][nH, z] = dG0
        else:
            self.comm.execute("CREATE TABLE gc_cid2prm (cid INT, nH INT, z INT, dG0 REAL);")
            for cid in self.kegg().get_all_cids_with_inchi():
                try:
                    pmap = self.estimate_pmap_keggcid(cid)
                    self.cid2pmap[cid] = pmap
                    for ((nH, z), dG0) in pmap.iteritems():
                        self.comm.execute("INSERT INTO gc_cid2prm VALUES(?,?,?,?)", (cid, nH, z, dG0))
                except (GroupDecompositionError, GroupMissingTrainDataError, KeyError, kegg.KeggParseException) as e:
                    sys.stderr.write("Warning, C%05d has an error: %s\n" % (cid, str(e)))
            self.comm.commit()
            
    def get_cid2dG0(self, pH=default_pH, I=default_I, T=default_T, most_abundant=False):
        cid2dG0 = {}
        for cid in self.cid2pmap.keys():
            pmap = self.cid2pmap[cid]
            if (len(pmap) == 0):
                raise ValueError("C%05d has an empty pmap" % cid)
            cid2dG0[cid] = self.pmap_to_dG0(pmap, pH, I, T, most_abundant)
        return cid2dG0           
        
    def train(self, obs_fname, use_dG0_format=False):
        if (use_dG0_format):
            self.read_training_data_dG0(obs_fname)
        else:
            self.read_training_data(obs_fname)
        self.group_contributions = self.linear_regression_train()

        sys.stderr.write("Storing the group contribution data in the database ... \n")

        self.comm.execute("DROP TABLE IF EXISTS contribution")
        self.comm.execute("CREATE TABLE contribution (gid INT, name TEXT, protons INT, charge INT, dG0_gr REAL)")
        for i in range(len(self.group_contributions)):
            j = int(self.nonzero_groups[i])
            (name, protons, charge) = self.all_groups[j]
            self.comm.execute("INSERT INTO contribution VALUES(?,?,?,?,?)", (j, name, protons, charge, self.group_contributions[i]))
            
        self.comm.execute("DROP TABLE IF EXISTS observation")
        self.comm.execute("CREATE TABLE observation (cid INT, name TEXT, protons INT, charge INT, dG0_f REAL, use_for TEXT)")

        for cid in self.cid2pmap_obs.keys():
            if (cid in self.cid_test_set):
                use_for = 'test'
            else:
                use_for = 'train'
            for ((nH, z), dG0) in self.cid2pmap_obs[cid].iteritems():
                self.comm.execute("INSERT INTO observation VALUES(?,?,?,?,?,?)", (cid, self.kegg().cid2name(cid), nH, z, dG0, use_for))
                
        self.comm.commit()
            
    def load_groups(self, group_fname=None):
        self.list_of_groups = [] # a map of inchi string to the location in the group_list
        if (group_fname != None):
            sys.stderr.write("Loading the list of groups from %s into the database ... " % group_fname)
            group_csv_file = csv.reader(open(group_fname, 'r'))
            group_csv_file.next()
        
            self.comm.execute("DROP TABLE IF EXISTS groups;")
            self.comm.execute("CREATE TABLE groups (gid INT, name TEXT, protons INT, charge INT, smarts TEXT, focal_atoms TEXT, remark TEXT)")
            gid = 0
            for row in group_csv_file:
                (group_name, protons, charge, smarts, focal_atom_set, remark) = row
                try: # this is just to make sure the Smarts definition are proper
                    pybel.Smarts(smarts)
                except IOError:
                    raise Exception("Cannot parse SMARTS from line %d: %s" % (group_csv_file.line_num, smarts))
                
                self.comm.execute("INSERT INTO groups VALUES(?,?,?,?,?,?,?)", (gid, group_name, int(protons), int(charge), smarts, focal_atom_set, remark))
                gid += 1
            sys.stderr.write("[DONE]\n")
            self.comm.commit()

        sys.stderr.write("Reading the list of groups from the database ... ")
        for row in self.comm.execute("SELECT * FROM groups"):
            (gid, group_name, protons, charge, smarts, focal_atom_set, remark) = row
            
            if (focal_atom_set != ""): # otherwise, consider all the atoms as focal atoms
                focal_atoms = set([int(i) for i in focal_atom_set.split('|')])
            else:
                focal_atoms = None

            self.list_of_groups.append((gid, group_name, protons, charge, str(smarts), focal_atoms))
        sys.stderr.write("[DONE]\n")

        mol = pybel.readstring('smi', 'C') # the specific compound is meaningless since we only want the group names
        (groups, unassigned_nodes) = self.decompose(mol)
        self.all_groups = [(group_name, protons, charge) for (group_name, protons, charge, node_sets) in groups] + [('origin', 0, 0)] # add the 'origin' group (i.e. for the 0-bias of the linear regression)
        self.all_group_names = ["%s [H%d %d]" % (group_name, protons, charge) for (group_name, protons, charge) in self.all_groups]
        self.all_group_protons = pylab.array([protons for (group_name, protons, charge) in self.all_groups])
        self.all_group_charges = pylab.array([charge for (group_name, protons, charge) in self.all_groups])
    
    def does_table_exist(self, table_name):
        for row in self.comm.execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def load_concentrations(self):
        self.media_list = []
        if (self.does_table_exist('compound_abundance')):
            for row in self.comm.execute("SELECT media FROM compound_abundance GROUP BY media"):
                self.media_list.append(row[0])
                
            self.cid2conc = {}
            for row in self.comm.execute("SELECT cid, media, concentration FROM compound_abundance"):
                (cid, media, conc) = row
                self.cid2conc[(cid, media)] = conc # in [M]

    def get_concentration(self, cid, c0=default_c0, media=None):
        if (cid == 1): # the concentration of water must always be 1
            return 1
        if (media == None):
            return c0 # Standard conditions = 1 [M]
        return self.cid2conc.get((cid, media), c0)

    def get_concentration_list(self, cid, c0=default_c0):
        """
            return a list of pairs of media names and concentrations of the provided CID
        """
        c_list = []
        for media in self.media_list:
            if ((cid, media) in self.cid2conc):
                c_list.append((media, self.cid2conc[(cid, media)]))
        return c_list
    
    def decompose(self, mol, ignore_protonations=False):
        """
            The flag 'ignore_protonations' should be used when decomposing a compound with lacing protonation
            representation (for example, the KEGG database doesn't posses this information).
            If this flag is set to True, it overrides the '(C)harge sensitive' flag in the groups file (i.e. - *PC)
        """
        unassigned_nodes = set(range(len(mol.atoms)))
        groups = []
        for (gid, group_name, protons, charge, smarts_str, focal_atoms) in self.list_of_groups:
    
            if (group_name[0:2] == "*P"): # phosphate chains require a special treatment
                if (group_name[2] == "I" or ignore_protonations): # (I)gnore charges
                    pchain_groups = find_pchains(mol, ignore_protonations=True)
                elif (group_name[2] == "C"): # (C)harge sensitive
                    pchain_groups = find_pchains(mol, ignore_protonations=False)
                else:
                    raise Exception("Unrecognized phosphate wildcard: " + group_name)
                for (group_key, group_nodesets) in pchain_groups:
                    (group_name, protons, charge) = group_key
                    current_groups = []
                    for focal_set in group_nodesets:
                        if (focal_set.issubset(unassigned_nodes)): # check that the focal-set doesn't override an assigned node
                            current_groups.append(focal_set)
                            unassigned_nodes = unassigned_nodes - focal_set
                    groups.append((group_name, protons, charge, current_groups))
            else:
                current_groups = []
                for nodes in find_smarts(smarts_str, mol):
                    try:
                        if (focal_atoms != None):
                            focal_set = set([nodes[i] for i in focal_atoms])
                        else:
                            focal_set = set(nodes)
                    except IndexError:
                        sys.stderr.write("Focal set for group %s is out of range: %s" % (group_name, str(focal_atoms)))
                        sys.exit(-1)

                    if (focal_set.issubset(unassigned_nodes)): # check that the focal-set doesn't override an assigned node
                        current_groups.append(focal_set)
                        unassigned_nodes = unassigned_nodes - focal_set
                groups.append((group_name, protons, charge, current_groups))
        
        for nodes in find_smarts("[H]", mol): # ignore the hydrogen atoms when checking which atom is unassigned
            unassigned_nodes = unassigned_nodes - set(nodes)
        
        return (groups, unassigned_nodes)
        
    def groups_to_string(self, groups):    
        group_strs = []
        for (group_name, protons, charge, node_sets) in groups:
            if (len(node_sets) > 0):
                group_strs.append('%s [H%d %d] x %d' % (group_name, protons, charge, len(node_sets)))
        return " | ".join(group_strs)
    
    def get_decomposition_str(self, mol):
        (groups, unassigned_nodes) = self.decompose(mol)
        if (len(unassigned_nodes) > 0):
            raise GroupDecompositionError("Unable to decompose %s into groups" % mol.title)
        else:
            self.groups_to_string(groups)
    
    def groups_to_vector(self, groups):
        return [len(node_sets) for (group_name, protons, charge, node_sets) in groups] + [1] # add 1 for the 'origin' group
    
    def get_groupvec(self, mol):
        (groups, unassigned_nodes) = self.decompose(mol)
        if (len(unassigned_nodes) > 0):
            raise GroupDecompositionError("Unable to decompose %s into groups" % mol.title)
        else:
            return self.groups_to_vector(groups)
        
    def groupvec2str(self, groupvec):
        group_strs = []
        for i in range(len(self.all_group_names)):
            if (groupvec[i] > 0):
                group_strs.append('%s x %d' % (self.all_group_names[i], groupvec[i]))
        return " | ".join(group_strs)
    
    def groupvec2charge(self, groupvec):
        return int(pylab.dot(groupvec, self.all_group_charges))

    def groupvec2protons(self, groupvec):
        return int(pylab.dot(groupvec, self.all_group_protons))            

    def get_protonated_groupvec(self, mol):
        (groups, unassigned_nodes) = self.decompose(mol, ignore_protonations=True)
        if (len(unassigned_nodes) > 0):
            raise GroupDecompositionError("Unable to decompose %s into groups" % mol.title)
        else:
            # 'group_name_to_index' is a map from each group name to its indices in the groupvec
            # note that some groups appear more than once (since they can have multiple protonation
            # levels).
            group_name_to_index = {}

            # 'group_name_to_count' is a map from each group name to its number of appearences in 'mol'
            group_name_to_count = {}
            for i in range(len(groups)):
                (group_name, protons, charge, node_sets) = groups[i]
                group_name_to_index[group_name] = group_name_to_index.get(group_name, []) + [i]
                group_name_to_count[group_name] = group_name_to_count.get(group_name, 0) + len(node_sets)
            
            index_vector = [] # maps the new indices to the original ones that are used in groupvec

            # a list of pairs, each containing the 'count' of each group and the number of possible protonations.
            total_slots_pairs = [] 

            for group_name in group_name_to_index.keys():
                groupvec_indices = group_name_to_index[group_name]
                index_vector += groupvec_indices
                total_slots_pairs.append((group_name_to_count[group_name], len(groupvec_indices)))

            # generate all possible assignments of protonations. Each group can appear several times, and we
            # can assign a different protonation level to each of the instances.
            groupvec_list = []
            for assignment in multi_distribute(total_slots_pairs):
                v = [0] * len(index_vector)
                for i in range(len(v)):
                    v[index_vector[i]] = assignment[i]
                groupvec_list.append(v + [1]) # add 1 for the 'origin' group
            return groupvec_list
            
    def get_pseudoisomers(self, mol):
        pseudoisomers = set()
        for groupvec in self.get_protonated_groupvec(mol):
            nH = self.groupvec2protons(groupvec)
            z = self.groupvec2charge(groupvec)
            pseudoisomers.add((nH, z))
        return sorted(list(pseudoisomers))
        
    def cid2pseudoisomers(self, cid):
        try:
            comp = self.kegg().cid2compound(cid)
            return self.get_pseudoisomers(comp.get_mol())
        except GroupDecompositionError:
            return [(self.kegg().cid2num_hydrogens(cid), self.kegg().cid2charge(cid))]
        except kegg.KeggParseException:
            return [(0, 0)]
    
    def mol2inchi(self, mol):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "inchi")
        return obConversion.WriteString(mol.OBMol).strip()        
    
    def read_training_data_dG0(self, obs_fname):
        """
            Finds all the compounds which have a valid dG0 in the dG0.csv file,
            and generates a regression matrix using the KEGG groups for these compounds.
            Return values is a tuple (X, y) where:
            X           - is the group regression matrix
            dG_obs      - is the observed dG0 vector.
        """
        X = []
        y = []
        self.mol_names = [] # is the list of Molecules (pybel class) used for the regression
        self.inchi2val_obs = {}
        self.cid2pmap_obs = {}
        self.cid_test_set = set()
        
        dG0_csv = csv.reader(open(obs_fname, 'r'))
        dG0_csv.next() # skip header row
    
        self.HTML.write('<h2><a name=compounds>List of compounds for training</a></h2>\n')
        self.HTML.write('Source File = %s<br>\n' % obs_fname)
        counter = 0
        for row in dG0_csv:
            (smiles, cid, compound_name, dG0, dH0, charge, hydrogens, Mg, use_for, ref, assumption) = row
            if (charge == ""):
                name = compound_name
            else:
                try:
                    name = "%s [%d]" % (compound_name, int(charge))
                except ValueError:
                    raise Exception("charge value is not an integer: " + charge)
            sys.stderr.write('Reading data for ' + name + '\n')
            self.HTML.write("<h3>%s, %s</h3>\n" % (name, ref))

            if (dG0 == ""):
                self.HTML.write('No data for &#x394;G<sub>f</sub><br>\n')
                continue

            if (use_for == "skip"):
                self.HTML.write('Compound marked as not to be used<br>\n')
                continue
                
            try:
                dG0 = float(dG0)
                self.HTML.write('&#x394;G<sub>f</sub> = %.2f<br>\n' % dG0)
            except ValueError:
                raise Exception("Invalid dG0: " + str(dG0))

            if (cid != ""):
                cid = int(cid)
                nH = int(hydrogens)
                z = int(charge)
                if (cid not in self.cid2pmap_obs):
                    self.cid2pmap_obs[cid] = {}
                self.cid2pmap_obs[cid][(nH, z)] = dG0

            if (use_for == "test"):
                self.cid_test_set.add(int(cid))
                self.HTML.write('Compound marked to be used only for testing (not training)<br>\n')
                continue
            elif (use_for == "train"):
                self.HTML.write('Compound marked to be used for training<br>\n')
            else:
                raise Exception("Unknown usage flag: " + use_for)

            if (smiles == ""):
                raise Exception("Cannot use compound '%s' for training if it lacks a SMILES string" % compound_name)
            try:
                self.HTML.write('SMILES = %s<br>\n' % smiles)
                mol = pybel.readstring('smiles', smiles)
            except TypeError:
                raise Exception("Invalid smiles: " + smiles)

            inchi = self.mol2inchi(mol)
            self.HTML.write('INCHI = %s<br>\n' % inchi)
            self.inchi2val_obs[inchi] = dG0
            mol.title = name
            
            img_fname = self.FIG_DIR + '/train_%05d.png' % counter
            counter += 1
            self.HTML.embed_img(img_fname, mol.title)
            self.HTML.write('<br>\n')
            try:
                mol.draw(show=False, filename=img_fname)
            except AssertionError:
                raise Exception("PyBel failed when trying to draw the compound %s" % compound_name)
    
            (groups, unassigned_nodes) = self.decompose(mol)
            if (len(unassigned_nodes) > 0):
                raise Exception("Compound '%s' could not be decomposed, I recommend moving it to the 'test only' set (SMILES = %s)" % (compound_name, smiles))
                
            groupvec = self.groups_to_vector(groups)
            self.mol_names.append(name)
            X.append(groupvec)
            y.append(dG0)
            self.HTML.write("Decomposition = %s<br>\n" % self.groups_to_string(groups))
            if (int(hydrogens) != self.groupvec2protons(groupvec)):
                self.HTML.write("ERROR: Hydrogen count doesn't match: explicit = %d, formula = %d<br>\n" % (int(hydrogens), self.groupvec2protons(groupvec)))
            if (int(charge) != self.groupvec2charge(groupvec)):
                self.HTML.write("ERROR: Charge doesn't match: explicit = %d, formula = %d<br>\n" % (int(charge), self.groupvec2charge(groupvec)))
                
        if (X == []):
            raise Exception("Could not use any of the groups in dG0.csv, aborting.")
        
        self.group_matrix = pylab.array(X)
        self.obs = pylab.array(y).T
        self.save_training_data()
    
    def read_training_data(self, obs_fname):
        X = []
        y = []
        self.mol_list = [] # is the list of Molecules (pybel class) used for the regression
        self.inchi2val_obs = {}
        self.cid2pmap_obs = {}
        
        self.HTML.write('<h2><a name=compounds>List of compounds for training</a></h2>\n')
        self.HTML.write('Source File = %s<br>\n' % obs_fname)
        Km_csv = csv.reader(open(obs_fname, 'r'))
        Km_csv.next()
        
        for row in Km_csv:
            (cid, obs_val) = row
            if (cid == "" or obs_val == ""):
                continue
            cid = int(cid)
            obs_val = float(obs_val)
            self.cid2pmap_obs[cid][(0, 0)] = obs_val
            try:
                mol = self.kegg().cid2mol(cid)
                self.HTML.write('<h3><a href="%s">C%05d - %s</a></h3>\n' % (self.kegg().cid2link(cid), cid, self.kegg().cid2name(cid)))
                self.HTML.write('Observed value = %f <br>\n' % obs_val)
                
                inchi = self.mol2inchi(mol)
                self.HTML.write('INCHI = %s<br>\n' % inchi)
                self.inchi2val_obs[inchi] = obs_val
                
                groupvec = self.get_groupvec(mol)
                self.mol_list.append(mol)
                X.append(groupvec)
                y.append(obs_val)
                self.HTML.write('Decomposition = %s <br>\n' % self.get_decomposition_str(mol))
            except GroupDecompositionError:
                self.HTML.write('Could not be decomposed<br>\n')
            except KeyError:
                self.HTML.write('Compound has no INCHI in KEGG<br>\n')

        self.group_matrix = pylab.array(X)
        self.obs = pylab.array(y).T
        self.save_training_data()

    def save_training_data(self):
        n_obs = self.group_matrix.shape[0]
        n_groups = self.group_matrix.shape[1]
        self.comm.execute("DROP TABLE IF EXISTS train_group_matrix")
        sql_command = "CREATE TABLE train_group_matrix (" + ",".join(["g%d REAL" % i for i in range(n_groups)]) + ")"
        self.comm.execute(sql_command)
        for j in range(n_obs):
            sql_command = "INSERT INTO train_group_matrix VALUES(" + ','.join(["?"]*n_groups) + ")"
            self.comm.execute(sql_command, self.group_matrix[j, :].tolist())
        
        self.comm.execute("DROP TABLE IF EXISTS train_observations")
        sql_command = "CREATE TABLE train_observations (obs REAL)"
        self.comm.execute(sql_command)
        for i in range(n_obs):
            self.comm.execute("INSERT INTO train_observations VALUES(%f)" % self.obs[i])
        
        self.comm.execute("DROP TABLE IF EXISTS train_groups")
        self.comm.execute("CREATE TABLE train_groups (name TEXT, protons INT, charge INT)")
        for (group_name, protons, charge) in self.all_groups:
            self.comm.execute("INSERT INTO train_groups VALUES(?,?,?)", (group_name, protons, charge))

        self.comm.execute("DROP TABLE IF EXISTS train_molecules")
        self.comm.execute("CREATE TABLE train_molecules (name TEXT)")
        for name in self.mol_names:
            self.comm.execute("INSERT INTO train_molecules VALUES(\'%s\')" % name)

        self.comm.commit()
            
    def load_training_data(self):
        X = []
        for row in self.comm.execute("SELECT * FROM train_group_matrix"):
            X.append(list(row))
        self.group_matrix = pylab.array(X)

        y = []
        for row in self.comm.execute("SELECT obs FROM train_observations"):
            y.append(row[0])
        self.obs = pylab.array(y)
        
        self.mol_names = []
        for row in self.comm.execute("SELECT name FROM train_molecules"):
            self.mol_names.append(row[0])

    def export_training_data(self, prefix):
        gmat_csv = csv.writer(open(prefix + "group_matrix.csv", "w"))
        gmat_csv.writerow(["compound name"] + self.all_group_names + ["observed dG0"])
        (n_comp, n_groups) = self.group_matrix.shape
        for i in range(n_comp):
            gmat_csv.writerow([self.mol_names[i]] + [x for x in self.group_matrix[i, :]] + [self.obs[i]])

        glist_csv = csv.writer(open(prefix + "groups.csv", "w"))
        glist_csv.writerow(("GROUP NAME", "PROTONS", "CHARGE"))
        for (group_name, protons, charge) in self.all_groups:
            glist_csv.writerow((group_name, protons, charge))
    
    def linear_regression_train(self):
        self.nonzero_groups = pylab.find(pylab.sum(self.group_matrix, 0) > 0)
        nonzero_group_mat = self.group_matrix[:, self.nonzero_groups]
        inv_corr_mat = pylab.pinv(pylab.dot(nonzero_group_mat.T, nonzero_group_mat))
        group_contributions = pylab.dot(pylab.dot(inv_corr_mat, nonzero_group_mat.T), self.obs)
        return group_contributions
    
    def groupvec2val(self, groupvec):
        if (self.group_contributions == None):
            raise Exception("You need to first Train the system before using it to estimate values")

        missing_groups = set(pylab.find(groupvec)).difference(set(self.nonzero_groups))
        if (len(missing_groups) == 0):
            groupvec = [groupvec[i] for i in self.nonzero_groups]
            return pylab.dot(groupvec, self.group_contributions)
        else:
            raise GroupMissingTrainDataError("can't estimate because some groups have no training data", missing_groups)
    
    def estimate_val(self, mol):
        try:
            groupvec = self.get_groupvec(mol)
        except GroupDecompositionError as e:
            inchi = self.mol2inchi(mol)
            if (inchi in self.inchi2val_obs):
                return self.inchi2val_obs[inchi]
            else:
                raise e
        
        return self.groupvec2val(groupvec)

    def write_regression_report(self, html_writer):
        group_matrix_reduced = self.group_matrix[:, self.nonzero_groups]
        html_writer.write('<h2><a name="regression">Regression</a></h2>\n')
        html_writer.write('<ul><li>%d compounds</li><li>%d groups</li><li>%d rank</li></ul>\n' % \
                        (group_matrix_reduced.shape[0], group_matrix_reduced.shape[1], matrixrank(group_matrix_reduced)))
        html_writer.write('<table border="1">\n<tr><td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td><td>Group Vector</td></tr>')
        for i in range(group_matrix_reduced.shape[0]):
            html_writer.write('<tr><td>%.2f</td><td>%s</td></tr>' % (self.obs[i], ' '.join([str(x) for x in group_matrix_reduced[i, :]])))
        html_writer.write('</table>')
        
        html_writer.write('<h2><a name="group_contrib">Group Contributions</a></h2>\n')
        html_writer.write('<table border="1">')
        html_writer.write('  <tr><td>Group Name</td><td>&#x394;<sub>gr</sub>G [kJ/mol]</td><td>Appears in compounds</td></tr>\n')
        j = 0
        for i in range(len(self.all_group_names)):
            if (i in self.nonzero_groups):
                contribution = self.group_contributions[j]
                compound_list_str = ' | '.join([self.mol_names[k] for k in pylab.find(group_matrix_reduced[:, j] > 0)])
                html_writer.write('  <tr><td>%s</td><td>%8.2f</td><td>%s</td></tr>\n' % (self.all_group_names[i], contribution, compound_list_str))
                j += 1
            else:
                html_writer.write('  <tr><td>%s</td><td>?</td><td></td></tr>\n' % (self.all_group_names[i]))
        html_writer.write('</table>\n')

    def analyze_training_set(self):
        self.write_regression_report(self.HTML)
        
        n_obs = len(self.obs)
        val_obs = []
        val_est = []
        val_err = []
        deviations = []

        nonzero_group_mat = self.group_matrix[:, self.nonzero_groups]
        compounds_with_unique_groups = [] 
        for i in range(n_obs):
            subset = pylab.array(range(0, i) + range((i + 1), n_obs))
            reduced_mat = nonzero_group_mat[subset, :]
            
            groups_in_i = set(pylab.find(nonzero_group_mat[i, :] > 0).tolist())
            groups_in_others = set(pylab.find(pylab.sum(reduced_mat, 0) > 0).tolist())
            if (not groups_in_i.issubset(groups_in_others)):
                compounds_with_unique_groups.append(i)
            else:
                inv_corr_mat = pylab.pinv(pylab.dot(reduced_mat.T, reduced_mat))
                group_contributions = pylab.dot(pylab.dot(inv_corr_mat, reduced_mat.T), self.obs[subset])
                
                estimation = pylab.dot(nonzero_group_mat[i, :], group_contributions)
                error = self.obs[i] - estimation
                val_obs.append(self.obs[i])
                val_est.append(estimation)
                val_err.append(error)
                deviations.append((abs(error), self.mol_names[i], self.obs[i], estimation, error))
        
        sys.stderr.write("Writing the table of estimation errors for each compound ... ")
        self.HTML.write('<h2><a name="error_table">Compound Estimation Error</a></h2>\n')
        self.HTML.write('<b>std(error) = %.2f kJ/mol</b>\n' % pylab.std(val_err))
        self.HTML.write('<table border="1">')
        self.HTML.write('  <tr><td>Compound Name</td><td>&#x394;<sub>f</sub>G<sub>obs</sub> [kJ/mol]</td><td>Error [kJ/mol]</td><td>Remark</td></tr>\n')
        deviations.sort(reverse=True)
        for (abs_err, mol_name, obs, est, err) in deviations:
            self.HTML.write('  <tr><td>%s</td><td>%8.2f</td><td>%8.2f</td><td>%s</td></tr>\n' % (mol_name, obs, err, ""))
        self.HTML.write('</table>\n')
        sys.stderr.write('[DONE]\n')
        
        sys.stderr.write("Plotting graphs for observed vs. estimated ... ")
        pylab.figure()
        pylab.plot(val_obs, val_est, '.')
        pylab.xlabel('Observed (obs)')
        pylab.ylabel('Estimated (est)')
        pylab.hold(True)
        for (abs_err, mol_name, obs, est, err) in deviations:
            pylab.text(obs, est, mol_name, fontsize=4)
        pylab.savefig('%s/obs_vs_est.svg' % self.FIG_DIR, format='svg')
        pylab.savefig('%s/obs_vs_est.pdf' % self.FIG_DIR, format='pdf')
        self.HTML.write('<h3><a name="obs_vs_est">Observed vs. Estimated</a></h3>\n')
        self.HTML.embed_svg('%s/obs_vs_est.svg' % self.EXP_NAME, name='obs_vs_est', width=1000, height=800)
        
        pylab.figure()
        pylab.plot(val_obs, val_err, '+')
        pylab.xlabel('Observed (obs)')
        pylab.ylabel('Estimation error (est - obs)')
        pylab.hold(True)
        for (abs_err, mol_name, obs, est, err) in deviations:
            pylab.text(obs, err, mol_name, fontsize=4)
    
        pylab.savefig('%s/obs_vs_err.svg' % self.FIG_DIR, format='svg', orientation='landscape')
        pylab.savefig('%s/obs_vs_err.pdf' % self.FIG_DIR, format='pdf', orientation='landscape')
        self.HTML.write('<h3><a name="obs_vs_err">Observed vs. Error</a></h3>\n')
        self.HTML.embed_svg('%s/obs_vs_err.svg' % self.EXP_NAME, name='obs_vs_err', width=1000, height=800)
        sys.stderr.write('[DONE]\n')
        self.HTML.flush()

    def kegg(self):
        if (self._kegg == None):
            self._kegg = kegg.Kegg(log_file=self.LOG_FILE)
        return self._kegg        

    def estimate_pmap(self, mol):
        try:
            all_groupvecs = self.get_protonated_groupvec(mol)
        except GroupDecompositionError as e:
            raise GroupDecompositionError(str(e) + " (%s)" % mol.title)

        all_missing_groups = []
        pmap = {}
        for groupvec in all_groupvecs:
            nH = self.groupvec2protons(groupvec)
            z = self.groupvec2charge(groupvec)
            try:
                dG0 = self.groupvec2val(groupvec)
                pmap[nH, z] = dG0
            except GroupMissingTrainDataError as e:
                s = "Species nH = %d, z = %d : " % (nH, z) + ", ".join([self.all_group_names[g] for g in e.missing_groups])
                all_missing_groups.append(s)
        
        if (len(pmap) == 0):
            raise GroupMissingTrainDataError("All species of %s have missing groups:" % mol.title, all_missing_groups)            

        return pmap

    def estimate_pmap_keggcid(self, cid):
        """
            returns a list of 3-tuples of dG0 (untransformed), nH and z.
            Each tuple represents one of the pseudoisomers.
        """

        if (cid == 80): # H+
            return {(0, 0):0}

        if (cid in self.cid2pmap_obs and (self.use_measured_train_values or cid in self.cid_test_set)):
            return self.cid2pmap_obs[cid]
        else:
            mol = self.kegg().cid2mol(cid)
            mol.title = "C%05d" % cid
            try:
                pmap = self.estimate_pmap(mol)
            except GroupMissingTrainDataError as e:
                if (cid in self.cid2pmap_obs):
                    return self.cid2pmap_obs[cid]
                else:
                    raise e

            return pmap
    
    @staticmethod
    def pmap_to_dG0(pmap, pH=default_pH, I=default_I, T=default_T, most_abundant=False):
        if (len(pmap) == 0):
            raise ValueError("Empty pmap given to 'pmap_to_dG0'")
        dG0_trans = pylab.array([Thermodynamics.transform(dG0, nH, z, pH, I, T) for ((nH, z), dG0) in pmap.iteritems()])
        if (most_abundant):
            return min(dG0_trans)
        else:
            return - R * T * log_sum_exp(-dG0_trans / (R * T))

    def estimate_pKa_keggcid(self, cid, charge, T=default_T):
        """
            Estimates the pKa of the compound.
            Ka is the equilibrium constant for the protonation reaction
            from the pseudoisomer with 'charge' to the pseudoisomer with 'charge'+1 
        """
        dG0_p0 = None
        dG0_p1 = None
        for (dG0, nH, z) in self.estimate_pmap_keggcid(cid):
            if (z == charge):
                dG0_p0 = dG0
            elif (z == charge + 1):
                dG0_p1 = dG0
        
        if (dG0_p0 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge))
        if (dG0_p1 == None):
            raise GroupMissingTrainDataError("cannot calculate dG0_f for C%05d pseudoisomer with charge %d" % (cid, charge + 1))
        
        return (dG0_p0 - dG0_p1) / (R * T * pylab.log(10))

    def estimate_dG0(self, mol, pH=default_pH, I=default_I, T=default_T): # T = temperature in K
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).            
        """
        pmap = self.estimate_pmap(mol)
        if (pH.__class__ == float and I.__class__ == float):
            return self.pmap_to_dG0(pmap, pH, I, T)      
        else:
            if (pH.__class__ == float):
                pH = [pH]
            if (I.__class__ == float):
                I = [I]
    
            dG0_matrix = pylab.zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_matrix[i, j] = self.pmap_to_dG0(pmap, pH[i], I[j], T)      
            
            return dG0_matrix
            
    def estimate_dG0_keggcid(self, cid, pH=default_pH, I=default_I, T=default_T, most_abundant=False): # T = temperature in K
        """
            Calculates the standard transformed Gibbs energy of formation of the pseudoisomer group
            (according to Alberty).            
        """
        pmap = self.estimate_pmap_keggcid(cid)
        if (pH.__class__ == float and I.__class__ == float):
            return self.pmap_to_dG0(pmap, pH, I, T, most_abundant)      
        else:
            if (pH.__class__ == float):
                pH = [pH]
            if (I.__class__ == float):
                I = [I]
    
            dG0_matrix = pylab.zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_matrix[i, j] = self.pmap_to_dG0(pmap, pH[i], I[j], T, most_abundant)      
            
            return dG0_matrix
    
    def estimate_dG0_reaction(self, sparse_reaction, pH=default_pH, I=default_I, T=default_T, most_abundant=False):
        """
            Calculate the dG0 (standard Gibbs free energy change)
            pH and I can be either floats or lists of values. If both are floats, returns a float.
            If either pH and/or I are lists, returns a matrix.
            T must be a float.
            
            dG = dG0 + R T \sum_i{s_i * ln([C_i])}
            
            s_i - stoichiometric coeff of compound i
            [C_i] - concentration of compound C_i
        """
        
        pmaps = [self.estimate_pmap_keggcid(cid) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()

        if (pH.__class__ == float and I.__class__ == float):
            dG0_vector = [self.pmap_to_dG0(pmap, pH, I, T, most_abundant) for pmap in pmaps]
            dG0 = pylab.dot(stoichiometry_vector, dG0_vector)
            return dG0
        else:
            if (pH.__class__ == float or pH.__class__ == int):
                pH = [float(pH)]
            if (I.__class__ == float or I.__class__ == int):
                I = [float(I)]
    
            dG0_matrix = pylab.zeros((len(pH), len(I)))
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0_vector = [self.pmap_to_dG0(pmap, pH[i], I[j], T, most_abundant) for pmap in pmaps]
                    dG0_matrix[i, j] = pylab.dot(stoichiometry_vector, dG0_vector)
            
            return dG0_matrix

    def estimate_dG_reaction(self, sparse_reaction, pH=default_pH, I=default_I, T=default_T, c0=default_c0, media=None, most_abundant=False):
        """
            standard = False means to use the known concentrations of reactants as well
        """
        
        dG0 = self.estimate_dG0_reaction(sparse_reaction, pH, I, T, most_abundant)
        concentration_vector = [pylab.log(self.get_concentration(cid, c0, media)) for cid in sparse_reaction.keys()]
        stoichiometry_vector = sparse_reaction.values()
        
        return dG0 + R * T * pylab.dot(stoichiometry_vector, concentration_vector)        

    def estimate_dG_keggrid(self, rid, pH=default_pH, I=default_I, T=default_T, c0=default_c0, media=None, most_abundant=False):
        """
            Returns the transformed Gibbs free energy change of a reaction according to its RID.
            Can set the pH, I and T to non-standard values.
            When media == None, it means we should use standard conditions (i.e. dG0).
            When media == 'glucose' (for example), it uses the concentrations measured for growth on glucose media.
        """
        sparse_reaction = self.kegg().rid2sparse_reaction(rid) 
        try:
            return self.estimate_dG_reaction(sparse_reaction, pH, I, T, c0, media, most_abundant)
        except KeyError as e:
            raise KeyError("R%05d contains a compound which cannot be used\n" % rid + str(e))

    def estimate_dG0_reaction_formula(self, formula, pH=default_pH, I=default_I, T=default_T, most_abundant=False):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.estimate_dG0_reaction(sparse_reaction, pH, I, T, most_abundant)

    def estimate_dG_reaction_formula(self, formula, pH=default_pH, I=default_I, T=default_T, media=None, most_abundant=False):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.estimate_dG_reaction(sparse_reaction, pH, I, T, media, most_abundant)
        
    def cid2groupvec(self, cid):
        try:
            return self.get_groupvec(self.kegg().cid2mol(cid))
        except GroupDecompositionError:
            raise GroupDecompositionError("Unable to decompose %s (C%05d) into groups" % (self.kegg().cid2name(cid), cid))
            
    def rid2groupvec(self, rid):
        sparse_reaction = self.kegg().rid2sparse_reaction(rid) 
        group_matrix = pylab.matrix([self.cid2groupvec(cid) for cid in sparse_reaction.keys()])
        stoichiometry_vector = pylab.matrix(sparse_reaction.values())
        total_groupvec = pylab.dot(stoichiometry_vector, group_matrix)
        return total_groupvec.tolist()[0]

    def analyze_all_kegg_compounds(self, pH=[default_pH], I=[default_I], T=default_T, most_abundant=False):
        self.comm.execute("DROP TABLE IF EXISTS dG0_f;")
        self.comm.execute("CREATE TABLE dG0_f (cid INT, pH REAL, I REAL, T REAL, dG0 REAL);")
        self.comm.execute("DROP INDEX IF EXISTS dG0_f_idx;")
        self.comm.execute("CREATE UNIQUE INDEX dG0_f_idx ON dG0_f (cid, pH, I, T);")

        self.HTML.write('<h2><a name=kegg_compounds>&#x394;G<sub>f</sub> of KEGG compounds:</a></h2>')
        for cid in self.cid2pmap.keys():
            self.HTML.write('<p>\n')
            self.HTML.write('<h3>C%05d <a href="%s">[KEGG]</a></h3>\n' % (cid, self.kegg().cid2link(cid)))
            self.HTML.write('Name: %s<br>\n' % self.kegg().cid2name(cid))
            self.HTML.write('Formula: %s<br>\n' % self.kegg().cid2formula(cid))
            
            pmap = self.cid2pmap[cid]
            for i in range(len(pH)):
                for j in range(len(I)):
                    dG0 = self.pmap_to_dG0(pmap, pH[i], I[j], T, most_abundant)
                    self.comm.execute("INSERT INTO dG0_f VALUES (?,?,?,?,?)", (cid, pH[i], I[j], T, dG0))
                    self.HTML.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>f</sub> = %.2f kJ/mol or %.2f kcal/mol<br>\n' % (pH[i], I[j], T, dG0, dG0 / 4.2))
            self.HTML.write('</p>\n')
        self.comm.commit()
    
    def analyze_all_kegg_reactions(self, pH=[default_pH], I=[default_I], T=default_T, most_abundant=False):
        self.comm.execute("DROP TABLE IF EXISTS dG0_r;")
        self.comm.execute("CREATE TABLE dG0_r (rid INT, pH REAL, I REAL, T REAL, dG0 REAL);")
        self.comm.execute("DROP INDEX IF EXISTS dG0_r_idx;")
        self.comm.execute("CREATE UNIQUE INDEX dG0_r_idx ON dG0_r (rid, pH, I, T);")

        self.HTML.write('<h2><a name=kegg_compounds>&#x394;G<sub>r</sub> of KEGG reactions:</a></h2>')
        for rid in self.kegg().get_all_rids():
            self.HTML.write('<p>\n')
            self.HTML.write('<h3>R%05d <a href="%s">[KEGG]</a></h3>\n' % (rid, self.kegg().rid2link(rid)))
            try:
                self.HTML.write('Definition: %s<br>\n' % self.kegg().rid2reaction(rid).definition)
                self.HTML.write('Equation:   %s<br>\n' % self.kegg().rid2reaction(rid).equation)
                dG0 = self.estimate_dG_keggrid(rid, pH=pH, I=I, T=T, media=None, most_abundant=most_abundant)
                for i in range(len(pH)):
                    for j in range(len(I)):
                        self.comm.execute("INSERT INTO dG0_r VALUES (?,?,?,?,?)", (rid, pH[i], I[j], T, dG0[i, j]))
                        self.HTML.write('Estimated (pH=%f, I=%f, T=%f) &#x394;G\'<sub>r</sub> = %.2f kJ/mol<br>\n' % (pH[i], I[j], T, dG0[i, j]))
            except GroupDecompositionError as e:
                self.HTML.write('Warning, cannot decompose one of the compounds: ' + str(e) + '<br>\n')
            except GroupMissingTrainDataError as e:
                self.HTML.write('Warning, cannot estimate: ' + str(e) + '<br>\n')
            except kegg.KeggParseException as e:
                self.HTML.write('Warning, cannot parse INCHI: ' + str(e) + '<br>\n')
            except KeyError as e:
                self.HTML.write('Warning, cannot locate this reaction in KEGG: ' + str(e) + '<br>\n')
            self.HTML.write('</p>\n')
        self.comm.commit()

    def write_cid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(["CID"] + self.all_group_names)
        for cid in self.kegg().get_all_cids_with_inchi():
            try:
                groupvec = self.cid2groupvec(cid)
                csv_file.writerow([cid] + groupvec)
            except GroupDecompositionError as e:
                print str(e)
                continue
            except kegg.KeggParseException as e:
                print str(e)
                continue
            
    def write_rid_group_matrix(self, fname):
        csv_file = csv.writer(open(fname, 'w'))
        csv_file.writerow(["RID"] + self.all_group_names)
        group_names = self.all_group_names
        
        groupstr_to_counter = {}
        for rid in self.kegg().get_all_rids():
            try:
                groupvec = self.rid2groupvec(rid)
                csv_file.writerow([rid] + groupvec)
                groupstr = ""
                for i in range(len(group_names)):
                    if (groupvec[i] != 0):
                        groupstr += "%s : %d, " % (group_names[i], groupvec[i])
                groupstr_to_counter[groupstr] = groupstr_to_counter.get(groupstr, 0) + 1
            except GroupDecompositionError as e:
                print "R%05d: Cannot decompose at least one of the compounds\n" % rid + str(e)
            except kegg.KeggParseException as e:
                print str(e)
            except KeyError as e:
                print "R%05d: Cannot locate the reaction or one of the compounds in kegg\n" % rid + str(e)
        
        gstr_hist = [(counter, groupstr) for (groupstr, counter) in groupstr_to_counter.iteritems()]
        gstr_hist.sort(reverse=True)
        f = open("../res/groupstr_hist.txt", "w")
        for (counter, groupstr) in gstr_hist:
            f.write("%5d ; " % counter + groupstr + "\n")
        f.close()
        
    def save_contributions(self, filename):
        sys.stderr.write("Saving the group contribution data to: %s\n" % filename)
        csv_output = csv.writer(open(filename, "w"))
        csv_output.writerow(("row type", "ID", "name", "protons", "charge", "dG"))
        for i in range(len(self.group_contributions)):
            j = self.nonzero_groups[i]
            (name, protons, charge) = self.all_groups[j]
            csv_output.writerow(("CONTRIBUTION", j, name, protons, charge, self.group_contributions[i]))
            
        for cid in self.cid2pmap_obs:
            for ((nH, z), dG0) in self.cid2pmap_obs[cid].iteritems():
                if (cid in self.cid_test_set):
                    use_for = 'test'
                else:
                    use_for = 'train'
                csv_output.writerow(("OBSERVATION", cid, self.kegg().cid2name(cid), nH, z, dG0, use_for))
            
    def load_contributions(self):
        sys.stderr.write("Loading the group contribution data from the database ... ")
        self.nonzero_groups = []
        self.group_contributions = []
        self.cid2pmap_obs = {}
        self.cid_test_set = set()
        
        for row in self.comm.execute("SELECT * FROM contribution"):
            (gid, name, protons, charge, dG0_gr) = row
            self.nonzero_groups.append(gid)
            self.group_contributions.append(dG0_gr)
        
        for row in self.comm.execute("SELECT * FROM observation"):
            (cid, name, nH, z, dG0_f, use_for) = row
            self.cid2pmap_obs.setdefault(cid, {})
            self.cid2pmap_obs[cid][nH, z] = dG0_f
            if (use_for == 'test'):
                self.cid_test_set.add(cid)
                
        sys.stderr.write('[DONE]\n')
        
    def read_compound_abundance(self, filename):
        self.comm.execute("DROP TABLE IF EXISTS compound_abundance")
        self.comm.execute("CREATE TABLE compound_abundance (cid INT, media TEXT, concentration REAL)")
        csv_in = csv.reader(open(filename, 'r'))
        csv_in.next() # skip title row
        for row in csv_in:
            (compound, cid, glu, glu_min, glu_max, gly, gly_min, gly_max, ace, ace_min, ace_max, reference, use) = row
            if (cid == ""):
                continue
            if (use != '1'):
                continue
            try:
                self.comm.execute("INSERT INTO compound_abundance VALUES (?,?,?)", (int(cid), "glucose", float(glu)))
            except ValueError:
                pass
            try:
                self.comm.execute("INSERT INTO compound_abundance VALUES (?,?,?)", (int(cid), "glycerol", float(gly)))
            except ValueError:
                pass
            try:
                self.comm.execute("INSERT INTO compound_abundance VALUES (?,?,?)", (int(cid), "acetate", float(ace)))
            except ValueError:
                pass
        self.comm.commit()
        self.load_concentrations()
    
    def write_table_to_html(self, html_writer, table_name):
        column_names = []
        for row in self.comm.execute("PRAGMA table_info(%s)" % table_name):
            (index, name, type, can_be_null, default_value, stam) = row
            column_names.append(name)
        
        html_writer.write('<p><h2><a name="%s">Table %s:</a></h2>\n' % (table_name, table_name))
        html_writer.write('<table border="1">\n')
        html_writer.write('  <tr>' + "".join([('<td><b>' + str(s) + '</b></td>') for s in column_names]) + '</tr>\n')
        for row in self.comm.execute("SELECT * FROM %s" % table_name):
            html_writer.write('  <tr>' + "".join([('<td>' + str(s) + '</td>') for s in row]) + '</tr>\n')
        html_writer.write('</table>\n')
        html_writer.write('</p>\n')

    def write_query_to_html(self, html_writer, query, title=None, column_names=None):
        if (title == None):
            title = query
        html_writer.write('<p><h2>%s</h2>\n' % (title))
        html_writer.write('<table border="1">\n')
        if (column_names != None):
            html_writer.write('  <tr>' + "".join([('<td><b>' + str(s) + '</b></td>') for s in column_names]) + '</tr>\n')
        for row in self.comm.execute(query):
            html_writer.write('  <tr>' + "".join([('<td>' + str(s) + '</td>') for s in row]) + '</tr>\n')
        html_writer.write('</table>\n')
        html_writer.write('</p>\n')

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if (len(tokens) == 0):
            return default_value
        if (len(tokens) > 1):
            raise Exception("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    def get_reactions(self, module_name, field_map, html_writer):
        """
            read the list of reactions from the command file
        """
        
        if ("MODULE" in field_map):
            mid_str = field_map["MODULE"]
            if (mid_str[0] == 'M'):
                mid = int(mid_str[1:])
            else:
                mid = int(mid_str)
            (S, rids, fluxes, cids) = self.kegg().get_module(mid)
            html_writer.write('<li>Module <a href=http://www.genome.jp/dbget-bin/www_bget?M%05d>M%05d</a></li>\n' % (mid, mid))       
        else:
            (S, rids, fluxes, cids) = self.kegg().parse_explicit_module(field_map)

        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        if ("MAP_CID" in field_map):
            for line in field_map["MAP_CID"].split('\t'):
                (cid_before, cid_after) = [int(cid[1:]) for cid in line.split(None, 1)]
                if (cid_before in cids):
                    cids[cids.index(cid_before)] = cid_after
        
        return (S, rids, fluxes, cids)

    def write_reactions_to_html(self, html_writer, S, rids, fluxes, cids, show_cids=True):
        
        def vector_to_hypertext(v, cids, show_cids=True):
            sparse_reaction = {}
            for c in range(len(v)):
                sparse_reaction[cids[c]] = v[c]
            return self.sparse_reaction_to_hypertext(sparse_reaction, show_cids=show_cids)
        
        html_writer.write("<li>Reactions:</br><ul>\n")
        
        for r in range(S.shape[0]):
            html_writer.write('<li><a href=' + self.kegg().rid2link(rids[r]) + '>R%05d' % rids[r] + '</a>')
            html_writer.write(' : ' + vector_to_hypertext(S[r, :].flat, cids, show_cids=show_cids))
            if (fluxes[r] != 1):
                html_writer.write(' (x%g)' % fluxes[r])
            html_writer.write('</li>\n')
        
        v_total = pylab.dot(pylab.matrix(fluxes), S).flat
        html_writer.write('<li><b>Total:</b>  ' + vector_to_hypertext(v_total, cids, show_cids=show_cids) + '</li>\n')
        html_writer.write("</ul></li>\n")
    
    def analyze_profile(self, key, field_map, html_writer):
        html_writer.write('<p>\n')

        # embed the graph first (though the SVG file itself will be created only later in the code
        svg_fname_graph = '%s/%s_graph.svg' % (self.FIG_DIR, key)
        svg_fname_profile = '%s/%s_profile.svg' % (self.FIG_DIR, key)
        html_writer.embed_svg(svg_fname_profile, name=key, width=800, height=600)
        html_writer.write('</br>\n')
        html_writer.write('<ul>\n')
        html_writer.write('<li>Conditions:</br><ol>\n')
        # read the list of conditions from the command file
        conditions = []
        for condition in field_map["CONDITIONS"].split('\t'):
            (media, pH, I, T, c0) = (None, default_pH, default_I, default_T, default_c0)
            media = re.findall("media=([a-zA-Z_]+)", condition)[0]
            if (media == 'None'):
                media = None
            pH = GroupContribution.get_float_parameter(condition, "pH", default_pH)
            I = GroupContribution.get_float_parameter(condition, "I", default_I)
            T = GroupContribution.get_float_parameter(condition, "T", default_T)
            c0 = GroupContribution.get_float_parameter(condition, "c0", default_c0)
            conditions.append((media, pH, I, T, c0))
            html_writer.write('<li>Conditions: media = %s, pH = %g, I = %g M, T = %g K, c0 = %g</li>\n' % (media, pH, I, T, c0))
        html_writer.write('</ol></li>\n')
        
        # read the list of methods for calculating the dG
        methods = []
        if (kegg.parse_bool_field(field_map, 'MILO', True)):
            methods.append('MILO')
        if (kegg.parse_bool_field(field_map, 'ABUNDANT', False)):
            methods.append('ABUNDANT')
        if (kegg.parse_bool_field(field_map, 'HATZI', False)):
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

        html_writer.write('<li><a href=%s>Graph visualization</a></li>' % svg_fname_graph)
        
        (S, rids, fluxes, cids) = self.get_reactions(key, field_map, html_writer)
        self.write_reactions_to_html(html_writer, S, rids, fluxes, cids, show_cids=False)
        html_writer.write('</ul>')

        (Nr, Nc) = S.shape

        # draw a graph representation of the pathway        
        Gdot = self.kegg().draw_pathway(S, rids, cids)
        Gdot.write(svg_fname_graph, prog='dot', format='svg')

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = {}
        dG0_r = {}
        dG_f = {}
        dG_r = {}
        for (method, media, pH, I, T, c0, plot_key) in params_list:
            dG0_f[plot_key] = pylab.zeros((Nc, 1))
            dG_f[plot_key] = pylab.zeros((Nc, 1))
            for c in range(Nc):
                if (method == "MILO"):
                    dG0_f[plot_key][c] = self.estimate_dG0_keggcid(cids[c], pH=pH, I=I, T=T)
                elif (method == "ABUNDANT"):
                    dG0_f[plot_key][c] = self.estimate_dG0_keggcid(cids[c], pH=pH, I=I, T=T, most_abundant=True)
                elif (method == "HATZI"):
                    dG0_f[plot_key][c] = self.hatzi_dG0_keggcid(cids[c], pH=pH, I=I, T=T)
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
        pylab.figure()
        pylab.hold(True)
        data = pylab.zeros((Nr + 1, len(legend)))
        for i in range(len(legend)):
            for r in range(1, Nr + 1):
                data[r, i] = sum(dG_r[legend[i]][:r, 0])
        pylab.plot(data)
        pylab.legend(legend, loc="lower left")

        for i in range(len(rids)):
            pylab.text(i + 0.5, pylab.mean(data[i:(i + 2), 0]), rids[i], fontsize=6, horizontalalignment='center', backgroundcolor='white')
        
        pylab.xlabel("Reaction no.")
        pylab.ylabel("dG [kJ/mol]")
        pylab.savefig(svg_fname_profile, format='svg')
        html_writer.write('</p>')

    @staticmethod
    def dG_to_str(dG):
        if (pylab.isnan(dG)):
            return "N/A"
        else:
            return "%.1f" % dG
    
    @staticmethod    
    def pC_to_range(pC, c_mid=1e-3, ratio=3.0):
        return (c_mid * 10 ** (-ratio/(ratio+1.0) * pC), c_mid * 10 ** (1.0/(ratio+1.0) * pC))

    @staticmethod
    def conc_to_str(conc):
        if (pylab.isnan(conc)):
            return "N/A"
        else:
            return "%.2g" % conc
    
    def analyze_slack(self, key, field_map, html_writer):
        html_writer.write('<p>\n')
      
        # embed the graph first (though the SVG file itself will be created only later in the code
        svg_fname_graph = '%s/%s_graph.svg' % (self.FIG_DIR, key)
        svg_fname_slack = '%s/%s_slack.svg' % (self.FIG_DIR, key)
        html_writer.embed_svg(svg_fname_slack, name=key, width=800, height=600)

        html_writer.write('</br>\n')
        html_writer.write('<ul>\n')
        html_writer.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc) 
        c_mid = kegg.parse_float_field(field_map, 'C_MID', 1e-3)
        (pH, I, T) = (default_pH, default_I, default_T)
        concentration_bounds = deepcopy(self.kegg().cid2bounds)
        if ("CONDITIONS" in field_map):
            pH = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            I = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            T = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                concentration_bounds[cid] = (conc, conc)
                html_writer.write(', [C%05d] = %g\n' % (cid, conc))
            html_writer.write('</li>\n')
        html_writer.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        html_writer.write('<li><a href=%s>Graph visualization</a></li>\n' % svg_fname_graph)

        (S, rids, fluxes, cids) = self.get_reactions(key, field_map, html_writer)
        self.write_reactions_to_html(html_writer, S, rids, fluxes, cids, show_cids=False)
        html_writer.write('</ul>\n')
        
        physiological_pC = kegg.parse_float_field(field_map, "PHYSIO", 4)

        # draw a graph representation of the pathway        
        (Nr, Nc) = S.shape
        Gdot = self.kegg().draw_pathway(S, rids, cids)
        Gdot.write(svg_fname_graph, prog='dot', format='svg')

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = pylab.zeros((Nc, 1))
        ind_nan = []
        html_writer.write('<table border="1">\n')
        html_writer.write('  <td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f [kJ/mol]", "nH", "z"))
        for c in range(Nc):
            cid = cids[c]
            name = self.kegg().cid2name(cid)
            try:
                pmap = self.estimate_pmap_keggcid(cid)
                dG0_f[c] = self.pmap_to_dG0(pmap, pH, I, T)
                for ((nH, z), dG0) in pmap.iteritems():
                    html_writer.write('<tr><td>%05d</td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td>\n' % (cid, name, dG0, nH, z))
            
            except (kegg.KeggParseException, GroupMissingTrainDataError):
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = pylab.nan
                ind_nan.append(c)
                html_writer.write('<tr><td>%05d</td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td>\n' % (cid, name))
        html_writer.write('</table>\n')
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        pC = pylab.arange(0, 20, 0.1)
        B_vec = pylab.zeros(len(pC))
        #label_vec = [""] * len(pC)
        #limiting_reactions = set()
        for i in xrange(len(pC)):
            c_range = GroupContribution.pC_to_range(pC[i], c_mid=c_mid)
            (dG_f, concentrations, B) = find_feasible_concentrations(S, dG0_f, c_range, bounds=bounds)
            dG_r = pylab.dot(S, dG_f)
            
            B_vec[i] = B
            #curr_limiting_reactions = set(pylab.find(abs(dG_r - B) < 1e-9)).difference(limiting_reactions)
            #label_vec[i] = ", ".join(["%d" % rids[r] for r in curr_limiting_reactions]) # all RIDs of reactions that have dG_r = B
            #limiting_reactions |= curr_limiting_reactions

        try:
            (dG_f, concentrations, pCr) = find_pCr(S, dG0_f, c_mid, bounds=bounds)
        except LinProgNoSolutionException:
            pCr = None
            
        try:
            c_range = GroupContribution.pC_to_range(physiological_pC, c_mid=c_mid)
            (dG_f, concentrations, B_physiological) = find_feasible_concentrations(S, dG0_f, c_range, bounds) 
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
        
        pylab.figure()
        pylab.plot(pC, B_vec, 'b')
        #for i in xrange(len(pC)):
        #    pylab.text(pC[i], B_vec[i], label_vec[i], fontsize=6, horizontalalignment='left', backgroundcolor='white')
        pylab.xlabel('pC')
        pylab.ylabel('slack [kJ/mol]')
        (ymin, ymax) = pylab.ylim()
        (xmin, xmax) = pylab.xlim()
#        pylab.broken_barh([(xmin, pCr), (pCr, xmax)], (ymin, 0), facecolors=('yellow', 'green'), alpha=0.3)
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
        pylab.savefig(svg_fname_slack, format='svg')

        # write a table of the compounds and their dG0_f
        html_writer.write('<table border="1">\n')
        html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f' [kJ/mol]"))
        for c in range(Nc):
            compound = self.kegg().cid2compound(cids[c])
            cid_str = '<a href="%s">C%05d</a>' % (compound.get_link(), compound.cid)
            html_writer.write('<tr><td>%s</td><td>%s</td><td>%s</td>\n' % (cid_str, compound.name, GroupContribution.dG_to_str(dG0_f[c, 0])))
        html_writer.write('</table><br>\n')
        
        # write a table of the reactions and their dG0_r
        html_writer.write('<table border="1">\n')
        html_writer.write('  <td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG RID", "Reaction", "flux"))
        for r in range(Nr):
            rid_str = '<a href="http://www.genome.jp/dbget-bin/www_bget?rn:R%05d">R%05d</a>' % (rids[r], rids[r])
            spr = {}
            for c in pylab.find(S[r, :]):
                spr[cids[c]] = S[r, c]
            reaction_str = self.kegg().sparse_to_hypertext(spr)
            html_writer.write('<tr><td>%s</td><td>%s</td><td>%g</td>\n' % (rid_str, reaction_str, fluxes[r]))
        html_writer.write('</table><br>\n')
        
        html_writer.write('</p>\n')

    def analyze_margin(self, key, field_map, html_writer):
        html_writer.write('<p>\n')
      
        # embed the graph first (though the SVG file itself will be created only later in the code
        svg_fname_graph = '%s/%s_graph.svg' % (self.FIG_DIR, key)
        svg_fname_profile = '%s/%s_profile.svg' % (self.FIG_DIR, key)
        html_writer.embed_svg(svg_fname_profile, name=key, width=800, height=600)

        html_writer.write('</br>\n')
        html_writer.write('<ul>\n')
        html_writer.write('<li>Conditions:</br><ol>\n')
        # c_mid the middle value of the margin: min(conc) < c_mid < max(conc) 
        c_mid = kegg.parse_float_field(field_map, 'C_MID', 1e-3)
        (pH, I, T) = (default_pH, default_I, default_T)
        concentration_bounds = deepcopy(self.kegg().cid2bounds)
        if ("CONDITIONS" in field_map):
            pH = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "pH", default_pH)
            I = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "I", default_I)
            T = GroupContribution.get_float_parameter(field_map["CONDITIONS"], "T", default_T)
            html_writer.write('<li>Conditions: pH = %g, I = %g M, T = %g K' % (pH, I, T))
            for tokens in re.findall("C([0-9]+)=([0-9\.e\+\-]+)", field_map["CONDITIONS"]):
                cid = float(tokens[0])
                conc = float(tokens[1])
                concentration_bounds[cid] = (conc, conc)
                html_writer.write(', [C%05d] = %g\n' % (cid, conc))
            html_writer.write('</li>\n')
        html_writer.write('</ol></li>')
                    
        # The method for how we are going to calculate the dG0
        method = kegg.parse_string_field(field_map, "METHOD", "MILO")
        c_range = kegg.parse_vfloat_field(field_map, "C_RANGE", [1e-6, 1e-2])
        
        html_writer.write('<li><a href=%s>Graph visualization</a></li>\n' % svg_fname_graph)

        (S, rids, fluxes, cids) = self.get_reactions(key, field_map, html_writer)
        self.write_reactions_to_html(html_writer, S, rids, fluxes, cids, show_cids=False)
        html_writer.write('</ul>\n')

        # draw a graph representation of the pathway        
        (Nr, Nc) = S.shape
        Gdot = self.kegg().draw_pathway(S, rids, cids)
        Gdot.write(svg_fname_graph, prog='dot', format='svg')

        # calculate the dG_f of each compound, and then use S to calculate dG_r
        dG0_f = pylab.zeros((Nc, 1))
        ind_nan = []
        html_writer.write('<table border="1">\n')
        html_writer.write('  <td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n' % ("KEGG CID", "Compound Name", "dG0_f [kJ/mol]", "nH", "z"))
        for c in range(Nc):
            cid = cids[c]
            name = self.kegg().cid2name(cid)
            try:
                if (method == "MILO"):
                    pmap = self.estimate_pmap_keggcid(cid)
                elif (method == "HATZI"):
                    dG0 = self.hatzi_cid2dG[cid]
                    z = self.hatzi_cid2charge[cid]
                    pmap = {(z, z): dG0}
                else:
                    raise Exception("Unknown dG evaluation method: " + method)
                
                dG0_f[c] = self.pmap_to_dG0(pmap, pH, I, T)
                for ((nH, z), dG0) in pmap.iteritems():
                    html_writer.write('<tr><td>%05d</td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td>\n' % (cid, name, dG0, nH, z))
            
            except (kegg.KeggParseException, GroupMissingTrainDataError):
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = pylab.nan
                ind_nan.append(c)
                html_writer.write('<tr><td>%05d</td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td>\n' % (cid, name))
        html_writer.write('</table>\n')
        bounds = [concentration_bounds.get(cid, (None, None)) for cid in cids]
        res = {}
        try:
            res['pCr'] = find_pCr(S, dG0_f, c_mid=c_mid, ratio=3, bounds=bounds)
            #res['PCR2'] = find_unfeasible_concentrations(S, dG0_f, c_range, c_mid=c_mid, bounds=bounds)
            res['SLACK'] = find_feasible_concentrations(S, dG0_f, c_range, bounds=bounds)
        except LinProgNoSolutionException:
            html_writer.write('<b>No feasible solution found, cannot calculate the Margin</b>')
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 10
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 5
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        profile_fig = pylab.figure()
        profile_fig.hold(True)

        pylab.title('Thermodynamic profile', figure=profile_fig)
        pylab.ylabel('cumulative dG [kJ/mol]', figure=profile_fig)
        pylab.xlabel('Reaction KEGG ID', figure=profile_fig)
        pylab.xticks(pylab.arange(1, Nr + 1), ['R%05d' % rids[i] for i in xrange(Nr)], fontproperties=FontProperties(size=8), rotation=30)
        if (len(ind_nan) == 0): # TODO: if this is not the case, find a way to draw the dG0_r anyway
            dG0_r = pylab.dot(S, dG0_f)
            cum_dG0_r = pylab.cumsum([0] + [dG0_r[i, 0] * fluxes[i] for i in range(Nr)])
            pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG0_r, figure=profile_fig, label='Standard [1M]')
        pylab.grid(True, figure=profile_fig)

        for optimization in res.keys():
            (dG_f, conc, score) = res[optimization]
            if (score == None):
                continue

            dG_r = pylab.dot(S, dG_f)
            cum_dG_r = pylab.cumsum([0] + [dG_r[i, 0] * fluxes[i] for i in range(Nr)])
            pylab.plot(pylab.arange(0.5, Nr + 1), cum_dG_r, figure=profile_fig, label='%s = %.1f' % (optimization, score))

        pylab.legend()
        profile_fig.savefig(svg_fname_profile, format='svg')

        for optimization in res.keys():
            (dG_f, conc, score) = res[optimization]
            if (score == None):
                continue
            conc[ind_nan] = c_mid # give all compounds with unknown dG0_f the middle concentration value

            conc_fig = pylab.figure()
            conc_fig.suptitle('Concentrations (%s = %.1f)' % (optimization, score))
            pylab.xscale('log', figure=conc_fig)
            pylab.ylabel('Compound KEGG ID', figure=conc_fig)
            pylab.xlabel('Concentration [M]', figure=conc_fig)
            pylab.yticks(range(Nc, 0, -1), ["C%05d" % cid for cid in cids], fontproperties=FontProperties(size=8))
            pylab.plot(conc, range(Nc, 0, -1), '*b', figure=conc_fig)

            # TODO: instead of compound numbers, write the compound KEGG IDs

            x_min = conc.min() / 10
            x_max = conc.max() * 10
            y_min = 0
            y_max = Nc + 1
            
            for c in range(Nc):
                pylab.text(conc[c, 0] * 1.1, Nc - c, self.kegg().cid2name(cids[c]), \
                           figure=conc_fig, fontsize=6, rotation=0)
                (b_low, b_up) = bounds[c]
                if (b_low == None):
                    b_low = x_min
                if (b_up == None):
                    b_up = x_max
                pylab.plot([b_low, b_up], [Nc - c, Nc - c], '-k', linewidth=0.4)

            if (optimization == 'pCr'):
                c_range_opt = GroupContribution.pC_to_range(score, c_mid=c_mid, ratio=3.0)
                pylab.axvspan(c_range_opt[0], c_range_opt[1], facecolor='g', alpha=0.3, figure=conc_fig)
            else:
                pylab.axvspan(c_range[0], c_range[1], facecolor='g', alpha=0.3, figure=conc_fig)
            pylab.axis([x_min, x_max, y_min, y_max], figure=conc_fig)
            svg_fname_concentrations = '%s/%s_conc_%s.svg' % (self.FIG_DIR, key, optimization)
            conc_fig.savefig(svg_fname_concentrations, format='svg')
            html_writer.embed_svg(svg_fname_concentrations, name=key, width=800, height=600)

    def analyze_contour(self, key, field_map, html_writer):
        svg_fname = '%s/%s_contour.svg' % (self.FIG_DIR, key)
        html_writer.embed_svg(svg_fname, name=key, width=800, height=600)

        pH_list = kegg.parse_vfloat_field(field_map, "PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I_list = kegg.parse_vfloat_field(field_map, "I", [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        T = kegg.parse_float_field(field_map, "T", default_T)
        most_abundant = kegg.parse_bool_field(field_map, "ABUNDANT", False)
        formula = kegg.parse_string_field(field_map, "REACTION")
        c0 = kegg.parse_float_field(field_map, "C0", 1.0)
        media = kegg.parse_string_field(field_map, "MEDIA", "None")
        if (media == "None"):
            media = None
            
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        dG_r = self.estimate_dG_reaction(sparse_reaction, pH_list, I_list, T, c0, media, most_abundant)
        pylab.figure()
        
        pH_meshlist, I_meshlist = pylab.meshgrid(pH_list, I_list)
        CS = pylab.contour(pH_meshlist.T, I_meshlist.T, dG_r)       
        pylab.clabel(CS, inline=1, fontsize=10)
        pylab.xlabel("pH")
        pylab.ylabel("Ionic Strength")
        pylab.savefig(svg_fname, format='svg')
        html_writer.write('<br>\n' + self.sparse_reaction_to_hypertext(sparse_reaction) + '<br>\n')
        html_writer.write('</p>')

    def analyze_protonation(self, key, field_map, html_writer):
        svg_fname = '%s/%s_protonation.svg' % (self.FIG_DIR, key)
        html_writer.embed_svg(svg_fname, name=key, width=800, height=600)
        pH_list = kegg.parse_vfloat_field(field_map, "PH", [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        I = kegg.parse_float_field(field_map, "I", default_I)
        T = kegg.parse_float_field(field_map, "T", default_T)
        cid = kegg.parse_string_field(field_map, "COMPOUND")
        cid = int(cid[1:])
        pmap = self.estimate_pmap_keggcid(cid)
        data = pylab.zeros((len(pmap), len(pH_list)))
        for j in range(len(pH_list)):
            pH = pH_list[j]
            dG0_array = pylab.matrix([-Thermodynamics.transform(dG0, nH, z, pH, I, T) / (R * T) for ((nH, z), dG0) in pmap.iteritems()])
            dG0_array = dG0_array - max(dG0_array)
            p_array = pylab.exp(dG0_array)
            p_array = p_array / sum(p_array)
            data[:, j] = p_array    
        
        pylab.figure()
        pylab.plot(pH_list, data.T)
        prop = pylab.matplotlib.font_manager.FontProperties(size=10)
        name = self.kegg().cid2name(cid)
        pylab.legend(['%s [%d]' % (name, z) for ((nH, z), dG0) in pmap.iteritems()], prop=prop)
        pylab.xlabel("pH")
        pylab.ylabel("Pseudoisomer proportion")
        pylab.savefig(svg_fname, format='svg')
        pylab.savefig('%s/%s_protonation.pdf' % (self.FIG_DIR, key), format='pdf')
        html_writer.write('<table border="1">\n')
        html_writer.write('  <tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % ('dG0_f', '# hydrogen', 'charge'))
        for ((nH, z), dG0) in pmap.iteritems():
            html_writer.write('  <tr><td>%.2f</td><td>%d</td><td>%d</td></tr>\n' % (dG0, nH, z))
        html_writer.write('</table>')
        html_writer.write('</p>')

    def analyze_pathway(self, filename, html_fname):
        html_writer = HtmlWriter(html_fname)
        html_writer.write("<h1>Pathway analysis using Group Contribution Method</h1>\n")
        self.kegg()
        entry2fields_map = kegg.parse_kegg_file(filename)
        
        for key in sorted(entry2fields_map.keys()):
            sys.stderr.write("Analyzing pathway: " + key + " ... ")
            field_map = entry2fields_map[key]

            if (field_map.get("SKIP", "FALSE") == "TRUE"):
                sys.stderr.write("skipping\n")
                continue
            try:
                html_writer.write("<b>%s - %s</b>" % (field_map["NAME"], field_map["TYPE"]))
                html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % (key))
                html_writer.write('<div id="%s" style="display:none">' % key)
                html_writer.write('<h2>%s - %s</h2>\n' % (field_map["NAME"], field_map["TYPE"]))
            except KeyError:
                raise Exception("Both the 'NAME' and 'TYPE' fields must be defined for each pathway")

            if (field_map["TYPE"] == "PROFILE"):     
                self.analyze_profile(key, field_map, html_writer)
            elif (field_map["TYPE"] == "SLACK"):     
                self.analyze_slack(key, field_map, html_writer)
            elif (field_map["TYPE"] == "MARGIN"):     
                self.analyze_margin(key, field_map, html_writer)               
            elif (field_map["TYPE"] == "CONTOUR"):
                self.analyze_contour(key, field_map, html_writer)
            elif (field_map["TYPE"] == "PROTONATION"):
                self.analyze_protonation(key, field_map, html_writer)
            else:
                raise Exception("Unknown analysis type: " + field_map["TYPE"])
            html_writer.write('</div><br>\n')
            sys.stderr.write("[DONE]\n")
            
        html_writer.write('<h4>Measured concentration table:</h4>\n')
        html_writer.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % ('__concentrations__'))
        html_writer.write('<div id="%s" style="display:none">' % '__concentrations__')
        self.write_query_to_html(html_writer, \
            "select cid, media, 1000*concentration from compound_abundance ORDER BY cid, media", \
            title="Abundance", column_names=["CID", "Media", "concentration [mM]"])
        html_writer.write('</div><br>\n')
        html_writer.close()

    def add_vectors(self, stoichiometric_vector1, cid_vector1, stoichiometric_vector2, cid_vector2):
        cid_vector = list(set(cid_vector1 + cid_vector2))
        stoichiometric_vector = pylab.zeros(len(cid_vector))
        for i in range(len(cid_vector)):
            try:
                i1 = cid_vector1.index(cid_vector[i])
                stoichiometric_vector[i] += stoichiometric_vector1[i1]
            except ValueError:
                pass
            try:
                i2 = cid_vector2.index(cid_vector[i])
                stoichiometric_vector[i] += stoichiometric_vector2[i2]
            except ValueError:
                pass
        
        nonzero_values = pylab.find(stoichiometric_vector != 0)
        cid_vector = [cid_vector[i] for i in nonzero_values]
        stoichiometric_vector = stoichiometric_vector[nonzero_values]
        
        return (stoichiometric_vector, cid_vector)

    def sparse_reaction_to_hypertext(self, sparse_reaction, show_cids=True):
        s_left = []
        s_right = []
        for (cid, count) in sparse_reaction.iteritems():
            url = self.kegg().cid2link(cid)
            name = self.kegg().cid2name(cid)
            if (show_cids):
                show_string = "C%05d" % cid
                title = name + ': ' + ", ".join(["%s-%.2e[M]" % (m, c) for (m, c) in self.get_concentration_list(cid)])
            else:
                show_string = name
                title = "C%05d: " % cid + ", ".join(["%s-%.2e[M]" % (m, c) for (m, c) in self.get_concentration_list(cid)])
            
            if (count > 0):
                if (count == 1):
                    s_right.append('<a href="%s" title="%s">%s</a>' % (url, title, show_string))
                else:
                    s_right.append('%g <a href="%s" title="%s">%s</a>' % (count, url, title, show_string))
            elif (count < 0):
                if (count == -1):
                    s_left.append('<a href="%s" title="%s">%s</a>' % (url, title, show_string))
                else:
                    s_left.append('%g <a href="%s" title="%s">%s</a>' % (-count, url, title, show_string))
        return ' + '.join(s_left) + ' -> ' + ' + '.join(s_right)
    
    def load_hatzimanikatis_rid_data(self, csv_fname=None):
        if (csv_fname != None):
            csv_reader = csv.reader(open(csv_fname, 'r'))
            csv_reader.next()
            self.comm.execute("DROP TABLE IF EXISTS hatzi_rid_to_dG")
            self.comm.execute("CREATE TABLE hatzi_rid_to_dG (rid INT, dG REAL, uncertainty REAL)")
            for row in csv_reader:
                (rid, dG, uncertanty, reason) = row
                rid = int(rid[1:])
                try:
                    dG = float(dG) * 4.18400 # turn calories to Joules
                    uncertainty = float(uncertanty)
                except ValueError:
                    continue
                self.comm.execute("INSERT INTO hatzi_rid_to_dG VALUES(?,?,?)", (rid, dG, uncertainty))
            self.comm.commit()
            
        self.hatzi_rid2dG = {}
        if (self.does_table_exist('hatzi_rid_to_dG')):
            for row in self.comm.execute("SELECT rid, dG FROM hatzi_rid_to_dG"):
                (rid, dG) = row
                self.hatzi_rid2dG[rid] = dG

    def load_hatzimanikatis_cid_data(self, csv_fname=None):
        if (csv_fname != None):
            csv_reader = csv.reader(open(csv_fname, 'r'))
            csv_reader.next()
            self.comm.execute("DROP TABLE IF EXISTS hatzi_cid_to_dG")
            self.comm.execute("CREATE TABLE hatzi_cid_to_dG (cid INT, dG REAL, uncertainty REAL, charge INT)")
            self.comm.execute("INSERT INTO hatzi_cid_to_dG VALUES(?,?,?,?)", (80, 0.0, 0.0, 0)) # this causes the program to ignore H+ in reactions
            for row in csv_reader:
                (cid, dG, uncertanty, charge, reason) = row
                cid = int(cid[1:])
                try:
                    dG = float(dG) * 4.18400 # turn calories to Joules
                    uncertainty = float(uncertanty)
                    charge = int(charge)
                except ValueError:
                    continue
                self.comm.execute("INSERT INTO hatzi_cid_to_dG VALUES(?,?,?,?)", (cid, dG, uncertainty, charge))
            self.comm.commit()
            
        self.hatzi_cid2dG = {}
        self.hatzi_cid2charge = {}
        if (self.does_table_exist('hatzi_cid_to_dG')):
            for row in self.comm.execute("SELECT cid, dG, charge FROM hatzi_cid_to_dG"):
                (cid, dG, charge) = row
                self.hatzi_cid2dG[cid] = dG
                self.hatzi_cid2charge[cid] = charge

    def hatzi_dG0_keggcid(self, cid, pH=default_pH, I=default_I, T=default_T):
        dG0 = self.hatzi_cid2dG[cid]
        z = self.hatzi_cid2charge[cid]
        nH = z
        return Thermodynamics.transform(dG0, nH, z, pH, I, T)

    def hatzi_dG0_reaction(self, stoichiometry_vector, cid_vector, pH=default_pH, I=default_I, T=default_T):
        try:
            dG0_vector = pylab.array([self.hatzi_dG0_keggcid(cid, pH, I, T) for cid in cid_vector])
            dG0_r = pylab.dot(stoichiometry_vector, dG0_vector)
            return dG0_r
        except KeyError:
            return self.estimate_dG0_reaction(stoichiometry_vector, cid_vector, pH=pH, I=I, T=T)

    def hatzi_dG_reaction(self, sparse_reaction, pH=default_pH, I=default_I, T=default_T, c0=default_c0, media=None):
        try:
            dG0_vector = pylab.array([self.hatzi_dG0_keggcid(cid, pH, I, T) for cid in sparse_reaction.keys()])
            concentration_vector = [pylab.log(self.get_concentration(cid, c0, media)) for cid in sparse_reaction.keys()]
            stoichiometry_vector = sparse_reaction.values()
            dG0_r = pylab.dot(stoichiometry_vector, dG0_vector)
            return dG0_r + R * T * pylab.dot(stoichiometry_vector, concentration_vector)
        except KeyError:
            return self.estimate_dG_reaction(sparse_reaction, pH=pH, I=I, T=T)            

    def hatzi_dG0_formula(self, formula, pH=default_pH, I=default_I, T=default_T):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.hatzi_dG0_reaction(sparse_reaction, pH, I, T)

    def hatzi_dG_formula(self, formula, pH=default_pH, I=default_I, T=default_T):
        sparse_reaction = self.kegg().formula_to_sparse(formula)
        return self.hatzi_dG_reaction(sparse_reaction, pH, I, T)
    
    def analyze_decomposition_cid(self, cid):
        return self.analyze_decomposition(self.kegg().cid2mol(cid))

    def analyze_decomposition(self, mol):
        (groups, unassigned_nodes) = self.decompose(mol)
        s = ""
        s += "%30s | %2s | %2s | %s\n" % ("group name", "nH", "z", "nodes")
        s += "-" * 50 + "\n"
        for (group_name, protons, charge, node_sets) in groups:
            for n_set in node_sets:
                s += "%30s | %2d | %2d | %s\n" % (group_name, protons, charge, ','.join([str(i) for i in n_set]))
        if (len(unassigned_nodes) > 0):
            s += "\nUnassigned nodes: \n"
            s += "%2s | %2s | %2s\n" % ('#', 'AN', 'z')
            s += "-" * 20 + "\n"
            for i in unassigned_nodes:
                a = mol.atoms[i]
                s += "%2d | %2d | %2d\n" % (i, a.atomicnum, a.formalcharge)
        return s

#################################################################################################################
#                                                   MAIN                                                        #
#################################################################################################################
    
if (__name__ == '__main__'):

    if (False): # Affinity of substrate group contribution
        G = GroupContribution("Km", "../data/groups_martin.csv")
        G.train("../data/cid_vs_log-km.csv", use_dG0_format=False)
    
        G.write_cid_group_matrix("../res/group_matrix.csv")

    if (False): # for Shira
        G = GroupContribution("dG0", "../data/groups_for_shira.csv")
        
        mol = G.kegg().cid2mol(58)
        (groups, unassigned_nodes) = G.decompose(mol)
        for i in unassigned_nodes:
            print mol.atoms[i].type
        
        G.write_rid_group_matrix("../res/rid_to_groups.csv")


