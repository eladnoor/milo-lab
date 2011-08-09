from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter, NullHtmlWriter
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.util import lsum
from pygibbs.thermodynamic_constants import default_I, default_pH,\
    default_T
from pygibbs.kegg import Kegg
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy
from pygibbs.nist_regression import NistRegression
import logging
from toolbox.molecule import Molecule
from pygibbs.groups_data import GroupsData
from pygibbs.group_decomposition import GroupDecomposer, GroupDecompositionError
from pygibbs.dissociation_constants import DissociationConstants

def test_dissociation_table(diss, group_decomposer, id, ignore_missing_smiles=False):
        if diss is None:
            logging.warning('%s: does not appear in the dissociation table' % id)
            return
        
        nH, nMg = diss.GetMostAbundantPseudoisomer(pH=default_pH, I=default_I, pMg=14, T=default_T) 
        if nMg != 0:
            logging.warning('%s: default species has nMg = %d' % (id, nMg))
            return
        smiles = diss.GetSmiles(nH=nH, nMg=0)
        if not smiles:
            if not ignore_missing_smiles:
                logging.warning('%s: no SMILES in the dissociation table for nH = %d' % (id, nH))
            return

        logging.debug('%s: nH = %d, smiles = %s' % (id, nH, smiles))
        mol = Molecule.FromSmiles(smiles)
    
        try:
            decomposition = group_decomposer.Decompose(mol, ignore_protonations=False, strict=True)
        except GroupDecompositionError:
            return
        
        groupvec = decomposition.AsVector()
        logging.debug("%s: decomposition = %s" % (id, groupvec))
        gc_nH = decomposition.Hydrogens()
        if nH != gc_nH:
            logging.warning('%s: nH doesn\'t match: explicit = %d, decomposition = %d' % (id, nH, gc_nH))

def nist_dissociation_test():
    """
        Verifies that all the compounds in NIST are covered by the dissociation table, including SMILES strings.
    """
    db = SqliteDatabase('../res/gibbs.sqlite')
    nist_regression = NistRegression(db, html_writer=NullHtmlWriter())
    dissociation = nist_regression.dissociation
    groups_data = GroupsData.FromDatabase(db)
    group_decomposer = GroupDecomposer(groups_data)
    kegg = Kegg.getInstance()
    
    nist = nist_regression.nist
    for cid in nist.GetAllCids():
        id = "C%05d (%s)" % (cid, kegg.cid2name(cid))
        if kegg.cid2compound(cid).get_atom_bag() is None:
            logging.debug('%s: has no explicit formula' % id)
        else:
            diss = dissociation.GetDissociationTable(cid, create_if_missing=False)
            test_dissociation_table(diss, group_decomposer, id, ignore_missing_smiles=False)

def dissociation_decomposition_test():
    """
        Verifies that the decomposition of the compounds in the dissociation table match the nH of each species.
    """
    db = SqliteDatabase('../res/gibbs.sqlite')
    dissociation = DissociationConstants.FromFile()
    groups_data = GroupsData.FromDatabase(db)
    group_decomposer = GroupDecomposer(groups_data)
    kegg = Kegg.getInstance()

    for cid in dissociation.GetAllCids():
        id = "C%05d (%s)" % (cid, kegg.cid2name(cid))
        if kegg.cid2compound(cid).get_atom_bag() is None:
            logging.debug('%s: has no explicit formula' % id)
        else:
            diss = dissociation.GetDissociationTable(cid, create_if_missing=False)
            test_dissociation_table(diss, group_decomposer, id, ignore_missing_smiles=True)
        
def compare_charges():
    #db_public = SqliteDatabase('../data/public_data.sqlite')
    db_gibbs = SqliteDatabase('../res/gibbs.sqlite')
    print "Writing Compare Charges report to ../res/groups_report.html"
    html_writer = HtmlWriter("../res/groups_report.html")
    kegg = Kegg.getInstance()
    
    #pH, I, pMg, T = default_pH, default_I, default_pMg, default_T
    pH, I, pMg, T = default_pH, 0, 14, default_T
    
    cid2error = {}
    for row_dict in db_gibbs.DictReader("gc_errors"):
        cid = int(row_dict['cid'])
        cid2error[cid] = row_dict['error']

    estimators = {}
    estimators['hatzi'] = Hatzi(use_pKa=False)
    estimators['milo'] = PsuedoisomerTableThermodynamics.FromDatabase(
        db_gibbs, 'gc_pseudoisomers', name='Milo Group Contribution')
    
    all_cids = set(lsum([e.get_all_cids() for e in estimators.values()]))
    dict_list = []
    for cid in all_cids:
        try:
            name = kegg.cid2name(cid)
            link = kegg.cid2compound(cid).get_link()
        except KeyError:
            name = "unknown"
            link = ""
        row_dict = {'cid':'<a href="%s">C%05d</a>' % (link, cid),
                    'name':name, 'error':cid2error.get(cid, None)}
        for key, est in estimators.iteritems():
            try:
                pmap = est.cid2PseudoisomerMap(cid)
                dG0, dG0_tag, nH, z, nMg = pmap.GetMostAbundantPseudoisomer(pH, I, pMg, T)
            except MissingCompoundFormationEnergy:
                dG0, dG0_tag, nH, z, nMg = "", "", "", "", ""
            row_dict['nH_' + key] = nH
            row_dict['charge_' + key] = z
            row_dict['nMg_' + key] = nMg
            row_dict['dG0_' + key] = dG0
            row_dict['dG0_tag_' + key] = dG0_tag
        dict_list.append(row_dict)
        
    html_writer.write_table(dict_list, headers=['cid', 'name', 'charge_hatzi', 'charge_milo', 'error'])
    html_writer.close()
    
def main():
    #nist_dissociation_test()
    #dissociation_decomposition_test()
    compare_charges()

if __name__ == '__main__':
    main()
