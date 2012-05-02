#!/usr/bin/python

from pygibbs.kegg import Kegg
from pygibbs.group_decomposition import GroupDecompositionError, GroupDecomposer
from toolbox.molecule import Molecule
from pygibbs.groups_data import GroupsData
from pygibbs.group_vector import GroupVector
from pygibbs.dissociation_constants import DissociationConstants

def GetMolInput(dissociation):
    mols = [] # a list of pairs of Molecule objects and stoichiometric coefficients 
    while mols == []:
        print 'KEGG ID or SMILES (or Enter to quit):',
        s_input = raw_input()
        if not s_input:
            return []
        elif s_input[0] == 'C':
            try:
                cid = int(s_input[1:])
                mols = [(dissociation.GetMostAbundantMol(cid, pH=7, I=0.25, pMg=14, T=298.15), 1)]
                print "Compound:", mols[0][0].ToInChI()
            except ValueError:
                print 'syntax error: KEGG compound ID is bad (%s), please try again' % s_input
        elif s_input[0] == 'R':
            try:
                rid = int(s_input[1:])
                reaction = Kegg.getInstance().rid2reaction(rid)
                print "Reaction:", str(reaction)
                for cid, coeff in reaction.iteritems():
                    mols += [(dissociation.GetMostAbundantMol(cid, pH=7, I=0.25, pMg=14, T=298.15), coeff)]
            except ValueError:
                print 'syntax error: KEGG reaction ID is bad (%s), please try again' % s_input
        else:
            try:
                mols = [(Molecule.FromSmiles(s_input), 1)]
                print "Compound:", mols[0][0].ToInChI()
            except Exception:
                print 'unable to parse SMILES string, please try again'
        
    return mols

def DecomposeInputString(group_decomposer, dissociation, ignore_protonations=False):
    mols = GetMolInput(dissociation)
    
    if len(mols) == 0:
        return False
    if len(mols) == 1:
        try:
            decomposition = group_decomposer.Decompose(mols[0][0], ignore_protonations=ignore_protonations, strict=True)
            if ignore_protonations:
                all_groupvecs = decomposition.PseudoisomerVectors()
            else:
                all_groupvecs = [decomposition.AsVector()]
            for groupvec in all_groupvecs:
                print ">>> nH = %2d, z = %2d, nMg = %d : %s" % \
                    (groupvec.Hydrogens(), groupvec.NetCharge(), groupvec.Magnesiums(), str(groupvec))
        except GroupDecompositionError as e:
            print "Cannot decompose compound to groups: " + str(e)
            print e.GetDebugTable()
        mols[0][0].Draw()
    else:
        total_gv = GroupVector(group_decomposer.groups_data)
        for mol, coeff in mols:
            try:
                decomposition = group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
            except GroupDecompositionError as e:
                print "Cannot decompose compound to groups: " + str(e)
                continue
            gv = decomposition.AsVector()
            total_gv += gv * coeff
        print ">>> nH = %2d, z = %2d, nMg = %d : %s" % \
            (total_gv.Hydrogens(), total_gv.NetCharge(), total_gv.Magnesiums(), str(total_gv))
    
    return True

def main():
    #db = SqliteDatabase("../res/gibbs.sqlite")
    dissociation = DissociationConstants.FromPublicDB()
    fname = "../data/thermodynamics/groups_species.csv"
    groups_data = GroupsData.FromGroupsFile(fname, transformed=False)
    group_decomposer = GroupDecomposer(groups_data)
    
    while DecomposeInputString(group_decomposer, dissociation):
        pass

if __name__ == '__main__':
    main()
    
