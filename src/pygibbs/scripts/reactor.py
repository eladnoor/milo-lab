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

REACTOR_PATH = "../../Reactor/Debug/Reactor"


def getAllCids (kegg):
    cids = kegg.get_all_cids()
    inchis = map(lambda cid: kegg.cid2inchi(cid), cids)
    inchi2cid = {}
    for i, inchi in enumerate(inchis):
        inchi2cid[inchi] = cids[i]
    return inchi2cid

def cid2smile (kegg,cid):
    return kegg.cid2smiles(cid)
    

def smile2cid(inchi2cid, smile):
    return "C%05d"%inchi2cid[HashedMolecule(smile).getInchi()]

def runReactor (smiles_file, smirks_file, output_file):
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
                    html_writer=HtmlWriter('../res/mog.html'),
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
    r = Reaction.FromFormula("2 C00288 => C00048")
    pl.find_path("MOG", r)


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
