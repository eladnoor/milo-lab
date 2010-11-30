#!/usr/bin/python

import csv
import logging
import pybel
import sys


class GroupDecompositionError(Exception):
    pass


class GroupsDataError(Exception):
    pass


class MalformedGroupDefinitionError(GroupsDataError):
    pass


class GroupDecomposition(object):
    """Class representing the group decomposition of a molecule."""
    
    def __init__(self, mol, groups, unassigned_nodes):
        self.mol = mol
        self.groups = groups
        self.unassigned_nodes = unassigned_nodes

    def ToTableString(self):
        """Returns the decomposition as a tabular string."""
        spacer = '-' * 50 + '\n'
        l = ['%30s | %2s | %2s | %s\n' % ("group name", "nH", "z", "nodes"),
             spacer]
                
        for group_name, protons, charge, node_sets in self.groups:
            for n_set in node_sets:
                s = '%30s | %2d | %2d | %s\n' % (group_name, protons, charge,
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


class GroupDecomposer(object):
    """Decomposes compounds into their constituent groups."""
    
    def __init__(self, groups_data):
        """Construct a GroupDecomposer.
        
        Args:
            groups_data: a list of per-group data.
        """
        self.groups = groups_data

    @staticmethod
    def FromGroupsFile(filename):
        """Factory that initializes a GroupDecomposer from a CSV file."""
        assert filename
        list_of_groups = []
        
        logging.info('Reading the list of groups from %s ... ' % filename)
        group_csv_file = csv.reader(open(filename, 'r'))
        group_csv_file.next() # Skip the header
    
        gid = 0
        for row in group_csv_file:
            try:
                (group_name, protons, charge, smarts, focal_atoms, remark) = row
                
                if focal_atoms:
                    focal_atoms = [int(c) for c in focal_atoms.split('|')]
                if protons:
                    protons = int(protons)
                if charge:
                    charge = int(charge)
                
                # Check that the smarts are good.
                pybel.Smarts(smarts)
                
                list_of_groups.append((gid, group_name, protons, charge, str(smarts), focal_atoms))
            except ValueError, msg:
                print msg
                raise GroupsDataError('Wrong number of columns (%d) in one of the rows in %s: %s' %
                                      (len(row), filename, str(row)))
            except IOError:
                raise GroupsDataError('Cannot parse SMARTS from line %d: %s' %
                                      (group_csv_file.line_num, smarts))
            
            gid += 1
        logging.info('Done reading groups data.')
        
        return GroupDecomposer(list_of_groups)

    @staticmethod
    def _IsPhosphate(group_name):
        return group_name.startswith('*P')
    
    @staticmethod
    def _IgnoreCharges(group_name):
        # (I)gnore charges
        return group_name[2] == 'I'
    
    @staticmethod
    def _ChargeSensitive(group_name):
        # (C)harge sensitive
        return group_name[2] == 'C'
    
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

    # Not intended to be mutable.
    DEFAULT_PHOSPHATE_GROUP_MAP = {
            # Initial group
            ("-OPO3-", 1, 0): [],
            ("-OPO3-", 0, -1): [],
            # Middle group
            ("-OPO2-", 1, 0): [],
            ("-OPO2-", 0, -1): [],
            # Final group
            #("-OPO3",  2,  0): [],
            ("-OPO3", 1, -1): [],
            ("-OPO3", 0, -2): []
            }

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
            
        """
        group_map = dict(GroupDecomposer.DEFAULT_PHOSPHATE_GROUP_MAP)
        v_charge = [a.formalcharge for a in mol.atoms]
        
        for length in xrange(1, max_length + 1):
            
            # Find internal phosphate chains (ones in the middle of the molecule).
            smarts_str = GroupDecomposer._InternalPChainSmarts(length)
            for pchain in GroupDecomposer.FindSmarts(mol, smarts_str):
                if ignore_protonations:
                    group_map[("-OPO3-", 0, -1)].append(set(pchain[1:6]))
                else:
                    charge = v_charge[pchain[4]]
                    protons = charge + 1
                    group_map[("-OPO3-", protons, charge)].append(set(pchain[1:6]))
                
                for j in xrange(length-1):
                    pg_start = 6 + j*4  # start of the phosphate group
                    if (ignore_protonations):
                        group_map[("-OPO2-", 0, -1)].append(
                            set(pchain[pg_start:(pg_start + 4)]))
                    else:
                        charge = v_charge[pchain[pg_start + 2]]
                        protons = charge + 1
                        group_map[("-OPO2-", protons, charge)].append(
                            set(pchain[pg_start:(pg_start + 4)]))
        
            # Find terminal phosphate chains.
            smarts_str = GroupDecomposer._TerminalPChainSmarts(length)
            for pchain in GroupDecomposer.FindSmarts(mol, smarts_str):
                if ignore_protonations:
                    group_map[("-OPO3", 1, -1)].append(set(pchain[0:5]))
                else:
                    charge = v_charge[pchain[0]] + v_charge[pchain[3]]
                    protons = charge + 2
                    group_key = ("-OPO3", protons, charge)
                    if group_key not in group_map:
                        logging.warning('This protonation (%d) level is not allowed for terminal phosphate groups.' % protons)
                        logging.warning('Assuming the level is actually 1 - i.e. the charge is (-1).')
                        group_map[("-OPO3", 1, -1)].append(set(pchain[0:5]))
                    else:
                        group_map[group_key].append(set(pchain[0:5]))
                        
                for j in xrange(length-1):
                    pg_start = 5 + j*4  # start of the phosphate group
                    if (ignore_protonations):
                        group_map[("-OPO2-", 0, -1)].append(
                            set(pchain[pg_start:(pg_start + 4)]))
                    else:
                        charge = v_charge[pchain[pg_start + 2]]
                        protons = charge + 1
                        group_map[("-OPO2-", protons, charge)].append(
                            set(pchain[pg_start:(pg_start + 4)]))
                    
        return sorted(list(group_map.iteritems()))

    def Decompose(self, mol, ignore_protonations=False, assert_result=False):
        """
        Decompose a molecule into groups.
        
        The flag 'ignore_protonations' should be used when decomposing a compound with lacing protonation
        representation (for example, the KEGG database doesn't posses this information). If this flag is
        set to True, it overrides the '(C)harge sensitive' flag in the groups file (i.e. - *PC)
        
        Args:
            mol: the molecule to decompose.
            ignore_protonations: whether to ignore protonation levels.
            assert_result: whether to 
        """
        unassigned_nodes = set(range(len(mol.atoms)))
        groups = []
        
        for (gid, group_name, protons, charge, smarts_str, focal_atoms) in self.groups:
            # Phosphate chains require a special treatment
            if self._IsPhosphate(group_name):
                pchain_groups = None
                if self._IgnoreCharges(group_name) or ignore_protonations:
                    pchain_groups = self.FindPhosphateChains(mol, ignore_protonations=True)
                elif self._ChargeSensitive(group_name):
                    pchain_groups = self.FindPhosphateChains(mol, ignore_protonations=False)
                else:
                    raise MalformedGroupDefinitionError(
                        'Unrecognized phosphate wildcard: %s' % group_name)
                    
                for group_key, group_nodesets in pchain_groups:
                    group_name, protons, charge = group_key
                    current_groups = []
                    
                    for focal_set in group_nodesets:
                        if focal_set.issubset(unassigned_nodes):
                            # Check that the focal-set doesn't override an assigned node
                            current_groups.append(focal_set)
                            unassigned_nodes = unassigned_nodes - focal_set
                    groups.append((group_name, protons, charge, current_groups))
                    
            else:  # Not a phosphate group
                current_groups = []
                for nodes in self.FindSmarts(mol, smarts_str):
                    try:
                        if focal_atoms:
                            focal_set = set([nodes[i] for i in focal_atoms])
                        else:
                            focal_set = set(nodes)
                    except IndexError:
                        logging.error('Focal set for group %s is out of range: %s'
                                      % (group_name, str(focal_atoms)))
                        sys.exit(-1)

                    if (focal_set.issubset(unassigned_nodes)): # check that the focal-set doesn't override an assigned node
                        current_groups.append(focal_set)
                        unassigned_nodes = unassigned_nodes - focal_set
                groups.append((group_name, protons, charge, current_groups))
        
        # Ignore the hydrogen atoms when checking which atom is unassigned
        for nodes in self.FindSmarts(mol, '[H]'): 
            unassigned_nodes = unassigned_nodes - set(nodes)
        
        decomposition = GroupDecomposition(mol, groups, unassigned_nodes)
        if assert_result and decomposition.unassigned_nodes:
            raise GroupDecompositionError('Unable to decompose %s into groups.\n%s' %
                                          (mol.title, decomposition.ToTableString()))
        
        return decomposition


if __name__ == '__main__':
    decomposer = GroupDecomposer.FromGroupsFile('../../data/thermodynamics/groups_species.csv')
    
    smiles = 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O'
    smiles = 'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O'
    #smiles = 'C(C1C(C(C(C(O1)O)O)O)O)O'
    mol = pybel.readstring('smiles', smiles)
    
    decomposition = decomposer.Decompose(mol)
    print decomposition.ToTableString()

    
