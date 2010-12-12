#!/usr/bin/python

import csv
import logging
import pybel
import pylab
import sys

from toolbox import util


class GroupDecompositionError(Exception):
    pass


class GroupsDataError(Exception):
    pass


class MalformedGroupDefinitionError(GroupsDataError):
    pass


class Group(object):
    """Representation of a single group."""
    
    def __init__(self, id, name, protons, charge, nMg,
                 smarts=None, focal_atoms=None):
        self.id = id
        self.name = name
        self.protons = protons
        self.charge = charge
        self.nMg = nMg
        self.smarts = smarts
        self.focal_atoms = focal_atoms

    def IsPhosphate(self):
        return self.name.startswith('*P')
    
    def IgnoreCharges(self):
        # (I)gnore charges
        return self.name[2] == 'I'
    
    def ChargeSensitive(self):
        # (C)harge sensitive
        return self.name[2] == 'C'
    
    def FocalSet(self, nodes):
        """Get the set of focal atoms from the match.
        
        Args:
            nodes: the nodes matching this group.
        
        Returns:
            A set of focal atoms.
        """
        if self.focal_atoms:
            return set([nodes[i] for i in self.focal_atoms])
        return set(nodes)
    
    def __str__(self):
        return '%s [H%d Z%d Mg%d]' % (self.name, self.protons, self.charge, self.nMg)
    
    def __eq__(self, other):
        """Enable == checking."""
        return (str(self.name) == str(other.name) and
                self.protons == other.protons and
                self.charge == other.charge and
                self.nMg == other.nMg)
    
    def __hash__(self):
        """We are HASHABLE!"""
        return hash((self.name, self.protons, self.charge, self.nMg))
    

class GroupsData(object):
    """Contains data about all groups."""
    
    ORIGIN = Group('Origin', 'Origin', protons=0, charge=0, nMg=0)
    
    # Phosphate groups need special treatment, so they are defined in code...
    # TODO(flamholz): Define them in the groups file.
    INITIAL_P_1 = Group('InitialP1', '-OPO3-', 1, 0, 0)
    INITIAL_P_2 = Group('InitialP1', '-OPO3-', 0, -1, 0)
    MIDDLE_P_1 = Group('MiddleP1', '-OPO2-', 1, 0, 0)
    MIDDLE_P_2 = Group('MiddleP2', '-OPO2-', 0, -1, 0)
    FINAL_P_1 = Group('FinalP1', '-OPO3', 1, -1, 0)
    FINAL_P_2 = Group('FinalP2', '-OPO3', 0, -2, 0)
    #FINAL_P_3 = ('-OPO3',  2,  0)
    INITIAL_2_PHOSPHATE = Group('Initial2PChain', '-OPO3-OPO2-', 0, -2, 0)
    MIDDLE_2_PHOSPHATE = Group('Middle2PChain', '-OPO2-OPO2-', 0, -2, 0)

    # Phosphate groups with magnesium ions.
    MIDDLE_2_P_MG = Group('Initial2PChain', '-OPO2-OPO2-', 0, 0, 1)
    FINAL_P_MG_2 = Group('Middle2PChainMg', '-OPO3', 0, 0, 1)
    
    DEFAULT_INTERIOR_P = MIDDLE_P_2  # Ignoring protonation, interior chain
    DEFAULT_EXTERIOR_P = FINAL_P_1   # Ignoring protonation, exterior chain
    
    DEFAULTS = {'-OPO3-': INITIAL_P_2,
                '-OPO2-': MIDDLE_P_2,
                '-OPO3': FINAL_P_1,
                '-OPO3-OPO2-': INITIAL_2_PHOSPHATE,
                '-OPO2-OPO2-': MIDDLE_2_PHOSPHATE}  
    
    PHOSPHATE_GROUPS = (
        INITIAL_P_1, INITIAL_P_2,
        MIDDLE_P_1, MIDDLE_P_2,
        INITIAL_2_PHOSPHATE,
        MIDDLE_2_PHOSPHATE,
        FINAL_P_1, FINAL_P_2,
        MIDDLE_2_P_MG, FINAL_P_MG_2)

    MIDDLE_PHOSPHATES_TO_MGS = ((MIDDLE_2_PHOSPHATE, MIDDLE_2_P_MG),)    
    FINAL_PHOSPHATES_TO_MGS = ((FINAL_P_2, FINAL_P_MG_2),)
    
    def __init__(self, groups, include_mg=True):
        """Construct GroupsData.
        
        Args:
            groups: a list of Group objects.
        """
        self.groups = groups
        self.all_groups = self._GetAllGroups(self.groups)
        self.all_group_names = [str(g) for g in self.all_groups]
        self.all_group_hydrogens = pylab.array([g.protons for g in self.all_groups])
        self.all_group_charges = pylab.array([g.charge for g in self.all_groups])
        self.all_group_mgs = pylab.array([g.nMg for g in self.all_groups])
    
    @staticmethod
    def _GetAllGroups(groups):
        all_groups = []
        
        for group in groups:
            # Expand phosphate groups.
            if group.IsPhosphate():
                all_groups.extend(GroupsData.PHOSPHATE_GROUPS)
            else:
                all_groups.append(group)
        
        # Add the origin.
        all_groups.append(GroupsData.ORIGIN)
        return all_groups
    
    @staticmethod
    def _ConvertFocalAtoms(focal_atoms_str):
        return [int(c) for c in focal_atoms_str.split('|')]
    
    @staticmethod
    def FromGroupsFile(filename):
        """Factory that initializes a GroupData from a CSV file."""
        assert filename
        list_of_groups = []
        
        logging.info('Reading the list of groups from %s ... ' % filename)
        group_csv_file = csv.reader(open(filename, 'r'))
        group_csv_file.next() # Skip the header
    
        gid = 0
        for row in group_csv_file:
            try:
                (group_name, protons, charge, smarts, focal_atoms, unused_remark) = row
                
                if focal_atoms:
                    focal_atoms = GroupsData._ConvertFocalAtoms(focal_atoms)
                if protons:
                    protons = int(protons)
                if charge:
                    charge = int(charge)
                
                # Check that the smarts are good.
                pybel.Smarts(smarts)
                
                mgs = 0
                group = Group(gid, group_name, protons, charge, mgs, str(smarts),
                              focal_atoms)
                list_of_groups.append(group)
            except ValueError, msg:
                print msg
                raise GroupsDataError('Wrong number of columns (%d) in one of the rows in %s: %s' %
                                      (len(row), filename, str(row)))
            except IOError:
                raise GroupsDataError('Cannot parse SMARTS from line %d: %s' %
                                      (group_csv_file.line_num, smarts))
            
            gid += 1
        logging.info('Done reading groups data.')
        
        return GroupsData(list_of_groups)    

    @staticmethod
    def FromDatabase(db):
        """Factory that initializes a GroupData from a DB connection."""
        logging.info('Reading the list of groups from the database.')
        
        list_of_groups = []
        for row in db.Execute('SELECT * FROM groups'):
            (gid, group_name, protons, charge, nMg, smarts, focal_atom_set, unused_remark) = row
            
            if focal_atom_set: # otherwise, consider all the atoms as focal atoms
                focal_atoms = GroupsData._ConvertFocalAtoms(focal_atom_set)
            else:
                focal_atoms = None

            list_of_groups.append(Group(gid, group_name, protons, charge, nMg, str(smarts), focal_atoms))
        logging.info('Done reading groups data.')
        
        return GroupsData(list_of_groups)
    
    def ToDatabase(self, db):
        """Write the GroupsData to the database."""
        logging.info('Writing GroupsData to the database.')
        
        db.CreateTable('groups', 'gid INT, name TEXT, protons INT, charge INT, nMg INT, smarts TEXT, focal_atoms TEXT, remark TEXT')
        for group in self.groups:
            focal_atom_str = '|'.join([str(fa) for fa in group.focal_atoms])
            db.Insert('groups', [group.id, group.name, int(group.protons), int(group.charge), 
                                 int(group.nMg), group.smarts, focal_atom_str, ''])

        logging.info('Done writing groups data into database.')

    def Index(self, gr):
        return self.all_groups.index(gr)
    
        
class GroupVector(list):
    """A vector of groups."""
    
    def __init__(self, groups_data, iterable=None):
        """Construct a vector.
        
        Args:
            groups_data: data about all the groups.
            iterable: data to load in to the vector.
        """
        self.groups_data = groups_data
        
        if iterable:
            self.extend(iterable)
    
    def __str__(self):
        """Return a string representation of this group vector."""
        group_strs = []
        
        for i, name in enumerate(self.groups_data.all_group_names):
            if self[i]:
                group_strs.append('%s x %d' % (name, self[i]))
        return " | ".join(group_strs)
    
    def __iadd__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] += other[i]

    def __isub__(self, other):
        for i in xrange(len(self.groups_data.all_group_names)):
            self[i] -= other[i]
            
    def __add__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result.append(self[i] + other[i])
        return result

    def __sub__(self, other):
        result = GroupVector(self.groups_data)
        for i in xrange(len(self.groups_data.all_group_names)):
            result.append(self[i] - other[i])
        return result
    
    def NetCharge(self):
        """Returns the net charge."""
        return int(pylab.dot(self, self.groups_data.all_group_charges))
    
    def Hydrogens(self):
        """Returns the number of protons."""
        return int(pylab.dot(self, self.groups_data.all_group_hydrogens))

    def Magnesiums(self):
        """Returns the number of Mg2+ ions."""
        return int(pylab.dot(self, self.groups_data.all_group_mgs))
    

class GroupDecomposition(object):
    """Class representing the group decomposition of a molecule."""
    
    def __init__(self, groups_data, mol, groups, unassigned_nodes):
        self.groups_data = groups_data
        self.mol = mol
        self.groups = groups
        self.unassigned_nodes = unassigned_nodes

    def ToTableString(self):
        """Returns the decomposition as a tabular string."""
        spacer = '-' * 50 + '\n'
        l = ['%30s | %2s | %2s | %3s | %s\n' % ("group name", "nH", "z", "nMg", "nodes"),
             spacer]
                
        for group, node_sets in self.groups:
            for n_set in node_sets:
                s = '%30s | %2d | %2d | %2d | %s\n' % (group.name, group.protons,
                                                       group.charge, group.nMg,
                                                       ','.join([str(i) for i in n_set]))
                l.append(s)

        if self.unassigned_nodes:
            l.append('\nUnassigned nodes: \n')
            l.append('%10s | %10s | %10s | %10s\n' % ('index', 'atomicnum',
                                                      'valence', 'charge'))
            l.append(spacer)
            
            for i in self.unassigned_nodes:
                a = self.mol.atoms[i]
                l.append('%10d | %10d | %10d | %10d\n' % (i, a.atomicnum,
                                                          a.heavyvalence, a.formalcharge))
        return ''.join(l)

    def __str__(self):
        """Convert the groups to a string."""        
        group_strs = []
        for group, node_sets in self.groups:
            if node_sets:
                group_strs.append('%s [H%d %d %d] x %d' % (group.name, group.protons,
                                                           group.charge, group.nMg,
                                                           len(node_sets)))
        return " | ".join(group_strs)
    
    def AsVector(self):
        """Return the group in vector format.
        
        Note: self.groups contains an entry for *all possible* groups, which is
        why this function returns consistent values for all compounds.
        """
        group_vec = GroupVector(self.groups_data)
        for unused_group, node_sets in self.groups:
            group_vec.append(len(node_sets))
        group_vec.append(1) # The origin
        return group_vec
    
    def NetCharge(self):
        """Returns the net charge."""
        return self.AsVector().NetCharge()
    
    def Hydrogens(self):
        """Returns the number of protons."""
        return self.AsVector().Hydrogens()
    
    def Magnesiums(self):
        """Returns the number of Mg2+ ions."""
        return self.AsVector().Magnesiums()
    
    def CountGroups(self):
        """Returns the total number of groups in the decomposition."""
        return sum([len(gdata[-1]) for gdata in self.groups])

    def PseudoisomerVectors(self):
        """Returns a list of group vectors, one per pseudo-isomer."""    
        if not self.CountGroups():
            logging.info('No groups in this decomposition, not calculating pseudoisomers.')
            return []
        
        # A map from each group name to its indices in the group vector.
        # Note that some groups appear more than once (since they can have
        # multiple protonation levels).
        group_name_to_index = {}

        # 'group_name_to_count' is a map from each group name to its number of appearances in 'mol'
        group_name_to_count = {}
        for i, gdata in enumerate(self.groups):
            group, node_sets = gdata
            group_name_to_index.setdefault(group.name, []).append(i)
            group_name_to_count[group.name] = group_name_to_count.get(group.name, 0) + len(node_sets)
        
        index_vector = [] # maps the new indices to the original ones that are used in groupvec

        # A list of per-group pairs (count, # possible protonation levels).
        total_slots_pairs = [] 

        for group_name, groupvec_indices in group_name_to_index.iteritems():
            index_vector += groupvec_indices
            total_slots_pairs.append((group_name_to_count[group_name],
                                      len(groupvec_indices)))

        # generate all possible assignments of protonations. Each group can appear several times, and we
        # can assign a different protonation level to each of the instances.
        groupvec_list = []
        for assignment in util.multi_distribute(total_slots_pairs):
            v = [0] * len(index_vector)
            for i in xrange(len(v)):
                v[index_vector[i]] = assignment[i]
            v += [1]  # add 1 for the 'origin' group
            groupvec_list.append(GroupVector(self.groups_data, v))
        return groupvec_list


class GroupDecomposer(object):
    """Decomposes compounds into their constituent groups."""
    
    def __init__(self, groups_data):
        """Construct a GroupDecomposer.
        
        Args:
            groups_data: a GroupsData object.
        """
        self.groups_data = groups_data

    @staticmethod
    def FromGroupsFile(filename):
        """Factory that initializes a GroupDecomposer from a CSV file."""
        assert filename
        gd = GroupsData.FromGroupsFile(filename)
        return GroupDecomposer(gd)
    
    @staticmethod
    def FromDatabase(db):
        """Factory that initializes a GroupDecomposer from a CSV file."""
        assert db
        gd = GroupsData.FromDatabase(db)
        return GroupDecomposer(gd)
    
    @staticmethod
    def FindSmarts(mol, smarts_str):
        """
        Corrects the pyBel version of Smarts.findall() which returns results as tuples,
        with 1-based indices even though Molecule.atoms is 0-based.

        Args:
            mol: the molecule to search in.
            smarts_str: the SMARTS query to search for.
        
        Returns:
            The re-mapped list of SMARTS matches.
        """
        shift_left = lambda m: [(n - 1) for n in m] 
        return map(shift_left, pybel.Smarts(smarts_str).findall(mol))

    @staticmethod
    def _InternalPChainSmarts(length):
        return ''.join(['CO', 'P(=O)([OH,O-])O' * length, 'C'])
    
    @staticmethod
    def _TerminalPChainSmarts(length):
        return ''.join(['[OH,O-]', 'P(=O)([OH,O-])O' * length, 'C'])

    @staticmethod
    def AttachMgToPhosphateChain(mol, chain_map, assigned_mgs):
        """Attaches Mg2+ ions the appropriate groups in the chain.
        
        Args:
            mol: the molecule.
            chain_map: the groups in the chain.
            assigned_mgs: the set of Mg2+ ions that are already assigned.
        
        Returns:
            The updated list of assigned Mg2+ ions. 
        """
        # For each Mg2+ we see, we attribute it to a phosphate group if
        # possible. We prefer to assign it to a terminal phosphate, but otherwise 
        # we assign it to a 'middle' group when there are 2 of them.
        def AddMg(p_group, pmg_group, mg):
            node_set = chain_map[p_group].pop(0)
            mg_index = mg[0]
            node_set.add(mg_index)
            assigned_mgs.add(mg_index)
            chain_map[pmg_group].append(node_set)
        
        all_pmg_groups = GroupsData.FINAL_PHOSPHATES_TO_MGS + GroupsData.MIDDLE_PHOSPHATES_TO_MGS
        for mg in GroupDecomposer.FindSmarts(mol, '[Mg+2]'):
            if mg[0] in assigned_mgs:
                continue
            
            for p_group, pmg_group in all_pmg_groups:
                if chain_map[p_group]:
                    AddMg(p_group, pmg_group, mg)
                    break

        return assigned_mgs

    @staticmethod
    def UpdateGroupMapFromChain(group_map, chain_map):
        """Updates the group_map by adding the chain."""
        for group, node_sets in chain_map.iteritems():
            group_map.get(group, []).extend(node_sets)
        return group_map

    @staticmethod
    def FindPhosphateChains(mol, max_length=3, ignore_protonations=False):
        """
        Chain end should be 'OC' for chains that do not really end, but link to carbons.
        Chain end should be '[O-1,OH]' for chains that end in an hydroxyl.
    
        Args:
            mol: the molecule to decompose.
            max_length: the maximum length of a phosphate chain to consider.
            ignore_protonations: whether or not to ignore protonation values.
        
        Returns:
            A list of 2-tuples (phosphate group, # occurrences).
        """
        group_map = dict((pg, []) for pg in GroupsData.PHOSPHATE_GROUPS)
        v_charge = [a.formalcharge for a in mol.atoms]
        assigned_mgs = set()
        
        def pop_phosphate(pchain, p_size):
            if len(pchain) < p_size:
                raise Exception('trying to pop more atoms than are left in the pchain')
            phosphate = pchain[0:p_size]
            charge = sum(v_charge[i] for i in phosphate)
            del pchain[0:p_size]
            return set(phosphate), charge
            
        def add_group(chain_map, group_name, charge, atoms):
            default = GroupsData.DEFAULTS[group_name]
            
            if ignore_protonations:
                chain_map[default].append(atoms)
            else:
                # NOTE(flamholz): We rely on the default number of magnesiums being 0 (which it is).
                protons = default.protons + charge - default.charge
                group = Group(default.id, group_name, protons, charge, default.nMg)
                if group not in chain_map:
                    logging.warning('This protonation (%d) level is not allowed for terminal phosphate groups.' % protons)
                    logging.warning('Using the default protonation level (%d) for this name ("%s").' %
                                    (default.protons, default.name))
                    chain_map[default].append(atoms)
                else:
                    chain_map[group].append(atoms)
        
        # For each allowed length
        for length in xrange(1, max_length + 1):
            # Find internal phosphate chains (ones in the middle of the molecule).
            smarts_str = GroupDecomposer._InternalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.iteritems())
            for pchain in GroupDecomposer.FindSmarts(mol, smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop() # Lose the carbons
                working_pchain.pop(0)
                
                if length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 5)
                    add_group(chain_map, '-OPO3-', charge, atoms)                    
                else:
                    atoms, charge = pop_phosphate(working_pchain, 9)
                    add_group(chain_map, '-OPO3-OPO2-', charge, atoms)
                
                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)
            
            assigned_mgs = GroupDecomposer.AttachMgToPhosphateChain(mol, chain_map,
                                                                    assigned_mgs)
            GroupDecomposer.UpdateGroupMapFromChain(group_map, chain_map)
            
            # Find terminal phosphate chains.
            smarts_str = GroupDecomposer._TerminalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.iteritems())
            for pchain in GroupDecomposer.FindSmarts(mol, smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop() # Lose the carbon
                
                atoms, charge = pop_phosphate(working_pchain, 5)
                add_group(chain_map, '-OPO3', charge, atoms)
                
                if not length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 4)
                    add_group(chain_map, '-OPO2-', charge, atoms)
                
                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)
                
            assigned_mgs = GroupDecomposer.AttachMgToPhosphateChain(mol, chain_map,
                                                                    assigned_mgs)
            GroupDecomposer.UpdateGroupMapFromChain(group_map, chain_map)

        return [(pg, group_map[pg]) for pg in GroupsData.PHOSPHATE_GROUPS]

    def Decompose(self, mol, ignore_protonations=False, strict=False):
        """
        Decompose a molecule into groups.
        
        The flag 'ignore_protonations' should be used when decomposing a compound with lacing protonation
        representation (for example, the KEGG database doesn't posses this information). If this flag is
        set to True, it overrides the '(C)harge sensitive' flag in the groups file (i.e. - *PC)
        
        Args:
            mol: the molecule to decompose.
            ignore_protonations: whether to ignore protonation levels.
            strict: whether to assert that there are no unassigned atoms.
        
        Returns:
            A GroupDecomposition object containing the decomposition.
        """
        unassigned_nodes = set(range(len(mol.atoms)))
        groups = []
        
        for group in self.groups_data.groups:
            # Phosphate chains require a special treatment
            if group.IsPhosphate():
                pchain_groups = None
                if group.IgnoreCharges() or ignore_protonations:
                    pchain_groups = self.FindPhosphateChains(mol, ignore_protonations=True)
                elif group.ChargeSensitive():
                    pchain_groups = self.FindPhosphateChains(mol, ignore_protonations=False)
                else:
                    raise MalformedGroupDefinitionError(
                        'Unrecognized phosphate wildcard: %s' % group.name)
                
                for phosphate_group, group_nodesets in pchain_groups:
                    current_groups = []
                    
                    for focal_set in group_nodesets:
                        if focal_set.issubset(unassigned_nodes):
                            # Check that the focal-set doesn't override an assigned node
                            current_groups.append(focal_set)
                            unassigned_nodes = unassigned_nodes - focal_set
                    groups.append((phosphate_group, current_groups))
                    
            else:  # Not a phosphate group
                current_groups = []
                for nodes in self.FindSmarts(mol, group.smarts):
                    try:
                        focal_set = group.FocalSet(nodes)
                    except IndexError:
                        logging.error('Focal set for group %s is out of range: %s'
                                      % (str(group), str(group.focal_atoms)))
                        sys.exit(-1)

                    # check that the focal-set doesn't override an assigned node
                    if focal_set.issubset(unassigned_nodes): 
                        current_groups.append(focal_set)
                        unassigned_nodes = unassigned_nodes - focal_set
                groups.append((group, current_groups))
        
        # Ignore the hydrogen atoms when checking which atom is unassigned
        for nodes in self.FindSmarts(mol, '[H]'): 
            unassigned_nodes = unassigned_nodes - set(nodes)
        
        decomposition = GroupDecomposition(self.groups_data, mol,
                                           groups, unassigned_nodes)
        
        if strict and decomposition.unassigned_nodes:
            raise GroupDecompositionError('Unable to decompose %s into groups.\n%s' %
                                          (mol.title, decomposition.ToTableString()))
        
        return decomposition


def main():
    decomposer = GroupDecomposer.FromGroupsFile('../data/thermodynamics/groups_species.csv')
    
    atp = 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O'
    coa = 'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O'
    glucose = 'C(C1C(C(C(C(O1)O)O)O)O)O'
    mgatp = 'C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-].[Mg+2].[Mg+2]'
    ctp = 'C1=CN(C(=O)N=C1N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O'

    smiless = [
               ('ATP', atp),
               ('CoA', coa), ('Glucose', glucose), ('MgAtp', mgatp),
               #('CTP', ctp)
               ]
    mols = [(name, pybel.readstring('smiles', s)) for name, s in smiless]

    for name, mol in mols:
        print name    
        decomposition = decomposer.Decompose(mol)
        print decomposition.ToTableString()
        print decomposition.AsVector()
        print 'Net charge', decomposition.NetCharge()
        print 'Hydrogens', decomposition.Hydrogens()
        
        for v in decomposition.PseudoisomerVectors():
            print v


if __name__ == '__main__':
    main()   
