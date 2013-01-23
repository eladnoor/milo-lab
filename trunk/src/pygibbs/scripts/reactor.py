###############################################################
#
# The script is run as follows:
# Given a file of list of smiles and a file of list of smirks,
# react over the smiles with each of the smirks.
# The output is then added to the pathfinder script run to find
# a path using all kegg reactions and the new reactions.
#
# Note! The script should run from the src directory.
# Make sure the the REACTOR_PATH variable is correct and point
# to the daylight compiled script.
#
###############################################################

from pygibbs.kegg import Kegg
import subprocess
from pathfinder.mol_list import HashedMolecule
from pygibbs.kegg_reaction import Reaction
from pygibbs.pathologic import Pathologic
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
from pygibbs.stoichiometric_lp import OptimizationMethods
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pathologic_examples import add_cofactor_reactions,add_redox_reactions

REACTOR_PATH = "./reactor/Debug/Reactor"


def getAllCids (kegg):
    "Get all the cids in kegg, given a Kegg object"
    cids = kegg.get_all_cids()
    inchis = map(lambda cid: kegg.cid2inchi(cid), cids)
    inchi2cid = {}
    for i, inchi in enumerate(inchis):
        inchi2cid[inchi] = cids[i]
    return inchi2cid

def cid2smile (kegg,cid):
    "Convert cid into smile, using Kegg object"
    return kegg.cid2smiles(cid)
    

def smile2cid(inchi2cid, smile):
    "Convert smile into cid through inchi, therefore using dictionary that converts inchi to cids"
    return "C%05d"%inchi2cid[HashedMolecule(smile).getInchi()]

def runReactor (smiles_file, smirks_file, output_file):
    """
    Run the daylight script that given a list of smiles as substrates and smirks as reactions,
    react on each substrate-reaction pair to get the product.
    smile_file  - filename with the list of all the molecules as smiles 
    smirks_file - filename with the list of all the reactions as smirks
    output_file - filename prefix to save the output to
    """
    return subprocess.call("%s -m %s -r %s -o %s"
                           % (REACTOR_PATH, smiles_file, smirks_file, output_file), shell=True) == 0
                           
def convertReactionIntoCid (inchi2cid,reaction):
    """
    Gets reaction as string and return it as cid.
    For example: CC>>C=C.O will convert into C80045 => C19503 + C00001
    In case of incorrect reaction, return empty string.
    """
    try:
        return " <=> ".join(map(lambda arr: " + ".join(map(lambda mol: smile2cid(inchi2cid,mol), arr.split('.'))), reaction.split('>>')))
    except:
        return ''
    
def runPathologic(thermo, reactionList):
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/mog_finder.html'),
                    thermo=thermo,
                    max_solutions=None,
                    max_reactions=15,
                    maximal_dG=-3.0,
                    thermodynamic_method=OptimizationMethods.GLOBAL,
                    update_file=None)
    add_cofactor_reactions(pl)
    add_redox_reactions(pl)
    for r in reactionList:
        pl.add_reaction(Reaction.FromFormula(r, "Auto generate #%s"%hash(r)))
    pl.delete_reaction(134)
    pl.delete_reaction(344)
    pl.delete_reaction(575)
    pl.delete_reaction(212)
    #pl.add_reaction(Reaction.FromFormula('C00149 + C00006 <=> C00036 + C00005 + C00080',
    #                                     'malate + NADP+ = oxaloacetate + NADPH',343))
    #pl.add_reaction(Reaction.FromFormula('C00222 + C00010 + C00006 <=> C00083 + C00005',
    #                                     'malonate-semialdehyde + CoA + NADP+ = malonyl-CoA + NADPH',740))
    r = Reaction.FromFormula("2 C00288 => C00048")
    pl.find_path("MOG_finder", r)
    

def runBeta2Alpha(thermo, reactionList):
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/Beta2Alpha.html'),
                    thermo=thermo,
                    max_solutions=None,
                    max_reactions=15,
                    maximal_dG=0.0,
                    thermodynamic_method=OptimizationMethods.GLOBAL,
                    update_file=None)
    add_cofactor_reactions(pl)
    add_redox_reactions(pl)
    for r in reactionList:
        pl.add_reaction(Reaction.FromFormula(r, "Auto generate #%s"%hash(r)))
    r = Reaction.FromFormula("C00099 => C01401")
    pl.find_path("Beta2Alpha", r)


if __name__ == '__main__':
    from sys import argv
    from time import time
    usage = "Usage: %s <smiles_file> <smirks_file>"%argv[0]
    if len(argv) != 3:
        print usage
        exit(1)
    
    kegg = Kegg.getInstance()
    inchi2cid = getAllCids(kegg)
    
    outputFile = "../res/out_%s"%time()
    if not runReactor(argv[1], argv[2], outputFile):
        print "Files are incorrect\n"+usage
        
    rxnFile = open(outputFile + '_reactions.txt')
    rxnList = set()
    for rxn in rxnFile:
        cid =convertReactionIntoCid(inchi2cid,rxn)
        if cid:
            rxnList.add(cid)
    rxnFile.close()
    
    estimators = LoadAllEstimators()
    thermo = estimators['UGC']
    thermo.SetConditions(pH=7.5, I=0.2)
    runPathologic(thermo,rxnList)
