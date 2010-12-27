#!/usr/bin/python

import csv
import logging
import pybel
import pylab

from toolbox import smarts_util


class GroupsDataError(Exception):
    pass


class MalformedGroupDefinitionError(GroupsDataError):
    pass


class _AllAtomsSet(object):
    """A set containing all the atoms: used for focal atoms sets."""
    
    def __contains__(self, elt):
        return True


class FocalSet(object):
    
    def __init__(self, focal_atoms_str):
        if not focal_atoms_str:
            raise ValueError(
                'You must supply a non-empty focal atom string.'
                ' You may use "None" or "All" in the obvious fashion.')
        
        self.str = focal_atoms_str
        self.focal_atoms_set = None
        prepped_str = self.str.strip().lower()
        
        if prepped_str == 'all':
            self.focal_atoms_set = _AllAtomsSet()
        elif prepped_str == 'none':
            self.focal_atoms_set = set()
        else:
            self.focal_atoms_set = set([int(c) for c in self.str.split('|')])
    
    def __str__(self):
        return self.str
    
    def __contains__(self, elt):
        return self.focal_atoms_set.__contains__(elt)


class Group(object):
    """Representation of a single group."""
    
    def __init__(self, id, name, hydrogens, charge, nMg,
                 smarts=None, focal_atoms=None):
        self.id = id
        self.name = name
        self.hydrogens = hydrogens
        self.charge = charge
        self.nMg = nMg
        self.smarts = smarts
        self.focal_atoms = focal_atoms

    def _IsHydrocarbonGroup(self):
        return self.name.startswith('*Hc')

    def _IsSugarGroup(self):
        return self.name.startswith('*Su')
    
    def _IsAromaticRingGroup(self):
        return self.name.startswith('*Ar')
    
    def _IsHeteroaromaticRingGroup(self):
        return self.name.startswith('*Har')

    def IsPhosphate(self):
        return self.name.startswith('*P')
    
    def IgnoreCharges(self):
        # (I)gnore charges
        return self.name[2] == 'I'
    
    def ChargeSensitive(self):
        # (C)harge sensitive
        return self.name[2] == 'C'
    
    def IsCodedCorrection(self):
        """Returns True if this is a correction for which hand-written code.
           must be executed.
        """
        return (self._IsHydrocarbonGroup() or
                self._IsAromaticRingGroup() or
                self._IsHeteroaromaticRingGroup())

    @staticmethod
    def _IsHydrocarbon(mol):
        """Tests if a molecule is a simple hydrocarbon."""
        if smarts_util.FindSmarts(mol, '[!C;!c]'):
            # If we find anything other than a carbon (w/ hydrogens)
            # then it's not a hydrocarbon.
            return 0
        return 1    

    @staticmethod
    def _CountAromaticRings(mol):
        expressions = ['c1cccc1', 'c1ccccc1']
        count = 0
        for smarts_str in expressions:
            count += len(smarts_util.FindSmarts(mol, smarts_str))
        return count
    
    @staticmethod
    def _CountHeteroaromaticRings(mol):
        expressions = ['a1aaaa1', 'a1aaaaa1']
        count = 0
        for smarts_str in expressions:
            for match in smarts_util.FindSmarts(mol, smarts_str):
                atoms = set([mol.atoms[i].atomicnum for i in match])
                atoms.discard(8)  # Ditch carbons
                if atoms:
                    count += 1
        return count

    def GetCorrection(self, mol):
        """Get the value of the correction for this molecule."""
        if self._IsHydrocarbonGroup():
            return self._IsHydrocarbon(mol)
        elif self._IsAromaticRingGroup():
            return self._CountAromaticRings(mol)
        elif self._IsHeteroaromaticRingGroup():
            return self._CountHeteroaromaticRings(mol)
        
        raise TypeError('This group is not a correction.')
    
    def FocalSet(self, nodes):
        """Get the set of focal atoms from the match.
        
        Args:
            nodes: the nodes matching this group.
        
        Returns:
            A set of focal atoms.
        """        
        focal_set = set()
        for i, node in enumerate(nodes):
            if i in self.focal_atoms:
                focal_set.add(node)            
        return focal_set
    
    def __str__(self):
        return '%s [H%d Z%d Mg%d]' % (self.name, self.hydrogens, self.charge, self.nMg)
    
    def __eq__(self, other):
        """Enable == checking.
        
        Only checks name, protons, charge, and nMg.
        """
        return (str(self.name) == str(other.name) and
                self.hydrogens == other.hydrogens and
                self.charge == other.charge and
                self.nMg == other.nMg)
    
    def __hash__(self):
        """We are HASHABLE!
        
        Note that the hash depends on the same attributes that are checked for equality.
        """
        return hash((self.name, self.hydrogens, self.charge, self.nMg))


class GroupsData(object):
    """Contains data about all groups."""
    
    ORIGIN = Group('Origin', 'Origin', hydrogens=0, charge=0, nMg=0)
    
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
        self.all_group_hydrogens = pylab.array([g.hydrogens for g in self.all_groups])
        self.all_group_charges = pylab.array([g.charge for g in self.all_groups])
        self.all_group_mgs = pylab.array([g.nMg for g in self.all_groups])
    
    def Count(self):
        return len(self.all_groups)
    
    count = property(Count)
    
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
        if not focal_atoms_str:
            return _AllAtomsSet()
        if focal_atoms_str.lower().strip() == 'none':
            return set()
        
        return set([int(c) for c in focal_atoms_str.split('|')])
    
    @staticmethod
    def FromGroupsFile(filename):
        """Factory that initializes a GroupData from a CSV file."""
        assert filename
        list_of_groups = []
        
        logging.info('Reading the list of groups from %s ... ' % filename)
        group_csv_file = csv.reader(open(filename, 'r'))
        group_csv_file.next() # Skip the header
    
        gid = 0
        for row in csv.DictReader(open(filename)):
            if row.get('SKIP', False):
                logging.warning('Skipping group %s', row.get('NAME'))
                continue
            
            try:
                group_name = row['NAME']
                protons = row['PROTONS']
                charge = row['CHARGE']
                smarts = row['SMARTS']
                focal_atoms = row['FOCAL_ATOMS']
                remark = row['REMARK']
                
                focal_atoms = FocalSet(focal_atoms)
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
            except KeyError, msg:
                logging.error(msg)
                raise GroupsDataError('Failed to parse row.')
            except ValueError, msg:
                logging.error(msg)
                raise GroupsDataError('Wrong number of columns (%d) in one of the rows in %s: %s' %
                                      (len(row), filename, str(row)))
            except IOError:
                raise GroupsDataError('Cannot parse SMARTS from line %d: %s' %
                                      (group_csv_file.line_num, smarts))
            
            gid += 1
        logging.info('Done reading groups data.')
        
        return GroupsData(list_of_groups)    

    @staticmethod
    def FromDatabase(db, filename=None):
        """Factory that initializes a GroupData from a DB connection.
        
        Args:
            db: a Database object.
            filename: an optional filename to load data from when
                it's not in the DB. Will write to DB if reading from file.
        
        Returns:
            An initialized GroupsData object.
        """
        logging.info('Reading the list of groups from the database.')
        
        if not db.DoesTableExist('groups') and filename:
            groups_data = GroupsData.FromGroupsFile(filename)
            groups_data.ToDatabase(db)
            return groups_data
        
        # Table should exist.
        list_of_groups = []
        for row in db.Execute('SELECT * FROM groups'):
            (gid, group_name, protons, charge, nMg, smarts, focal_atom_set, unused_remark) = row
            try:
                focal_atoms = FocalSet(focal_atom_set)
            except ValueError as e:
                raise ValueError('Group #%d (%s): %s' % (gid, group_name, str(e)))
            list_of_groups.append(Group(gid, group_name, protons, charge, nMg, str(smarts), focal_atoms))
        logging.info('Done reading groups data.')
        
        return GroupsData(list_of_groups)
    
    def ToDatabase(self, db):
        """Write the GroupsData to the database."""
        logging.info('Writing GroupsData to the database.')
        
        db.CreateTable('groups', 'gid INT, name TEXT, protons INT, charge INT, nMg INT, smarts TEXT, focal_atoms TEXT, remark TEXT')
        for group in self.groups:
            focal_atom_str = str(group.focal_atoms)
            db.Insert('groups', [group.id, group.name, int(group.hydrogens), int(group.charge), 
                                 int(group.nMg), group.smarts, focal_atom_str, ''])

        logging.info('Done writing groups data into database.')

    def Index(self, gr):
        try:
            return self.all_groups.index(gr)
        except ValueError:
            raise ValueError('group %s is not defined' % str(gr))