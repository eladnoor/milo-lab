#!/usr/bin/python

import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from toolbox.molecule import Molecule
from pygibbs.groups_data import GroupsData

def GetMolInput(kegg):
    mol = None
    while mol is None:
        print 'KEGG compound ID or SMILES:',
        s_input = raw_input()
        try:
            cid = int(s_input)
            mol = kegg.cid2mol(cid)
        except ValueError:
            try:
                mol = Molecule.FromSmiles(s_input)
            except Exception:
                print 'unable to parse SMILES string, please try again'
    print "InChI:", mol.ToInChI()
    return mol

def DecomposeInputString(group_decomposer):
    kegg = Kegg.getInstance()
    mol = GetMolInput(kegg)
    try:
        decomposition = group_decomposer.Decompose(mol, ignore_protonations=True, strict=True)
        all_groupvecs = decomposition.PseudoisomerVectors()
        for groupvec in all_groupvecs:
            print ">>> nH = %2d, z = %2d, nMg = %d : %s" % \
                (groupvec.Hydrogens(), groupvec.NetCharge(), groupvec.Magnesiums(), str(groupvec))
    except GroupDecompositionError as e:
        print "Cannot decompose compound to groups: " + str(e)
        print e.GetDebugTable()

    mol.Draw()

def main():
    #db = SqliteDatabase("../res/gibbs.sqlite")
    fname = "../data/thermodynamics/groups_species.csv"
    groups_data = GroupsData.FromGroupsFile(fname, transformed=False)
    group_decomposer = GroupDecomposer(groups_data)
    
    while True:
        DecomposeInputString(group_decomposer)        

if __name__ == '__main__':
    main()
    
