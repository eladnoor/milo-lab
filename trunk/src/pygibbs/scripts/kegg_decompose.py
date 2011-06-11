#!/usr/bin/python

import sys

from pygibbs import flags
from pygibbs.thermodynamic_constants import default_T
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from pygibbs.group_decomposition import GroupDecompositionError
from toolbox.molecule import Molecule

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
    return mol

def DecomposeInputString(G, kegg):
    mol = GetMolInput(kegg)
    try:
        decomposition = G.Mol2Decomposition(mol, ignore_protonations=True)
        all_groupvecs = decomposition.PseudoisomerVectors()
        for groupvec in all_groupvecs:
            print groupvec
            print groupvec.Hydrogens(), groupvec.NetCharge(), groupvec.Magnesiums()
        #print decomposition.ToTableString()
        #print 'nH =', decomposition.Hydrogens()
        #print 'z =', decomposition.NetCharge()
        #print 'nMg = ', decomposition.Magnesiums()
        pmap = G.Mol2PseudoisomerMap(mol, ignore_protonations=True)
        print pmap
        dG0, dG0_tag, nH, z, nMg = pmap.GetMostAbundantPseudoisomer(G.pH, G.I, G.pMg, G.T)
        print dG0, dG0_tag, nH, z, nMg
    except GroupDecompositionError as e:
        print "Cannot decompose compound to groups: " + str(e)

    mol.Draw()

def main():
    options, _ = flags.MakeOpts().parse_args(sys.argv)
    kegg = Kegg.getInstance()
    
    db = SqliteDatabase("../res/gibbs.sqlite")
    G = GroupContribution(db=db)
    G.init()
    G.c_mid = options.c_mid
    G.pH = options.ph
    G.pMg = options.pmg
    G.I = options.i_s
    G.T = default_T

    print ('Parameters: T=%f K, pH=%.2g, pMg=%.2g, '
           'I=%.2gmM, Median concentration=%.2gM' % (G.T, G.pH, G.pMg, G.I, G.c_mid))
    
    while True:
        DecomposeInputString(G, kegg)        

if __name__ == '__main__':
    main()
    
