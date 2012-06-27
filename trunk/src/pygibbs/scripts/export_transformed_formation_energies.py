from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
import logging
from pygibbs.kegg import Kegg
import csv

def main():
    pH, I, pMg, T = 7.0, 0.25, 14.0, 298.15
    
    dissociation = DissociationConstants.FromPublicDB()
    kegg = Kegg.getInstance()
    obs_fname = "../data/thermodynamics/formation_energies.csv"
    res_fname = "../res/formation_energies_transformed.csv"
    
    train_species = PsuedoisomerTableThermodynamics.FromCsvFile(obs_fname, label='testing')
    csv_out = csv.writer(open(res_fname, 'w'))
    csv_out.writerow(['cid','name',"dG'0",'pH','I','pMg','T','anchor','compound_ref','remark'])
    for cid in train_species.get_all_cids():
        pmap = train_species.cid2PseudoisomerMap(cid)
        source = train_species.cid2source_string[cid]
        pmatrix = pmap.ToMatrix() # ToMatrix returns tuples of (nH, z, nMg, dG0)
        if len(pmatrix) != 1:
            raise Exception("multiple training species for C%05d" % cid)
        nH, charge, nMg, dG0 = pmatrix[0]
        name = "%s (%d)" % (kegg.cid2name(cid), nH)
        logging.info('Adding the formation energy of %s', name)
        diss_table = dissociation.GetDissociationTable(cid, create_if_missing=True)
        if diss_table is None:
            raise Exception("%s [C%05d, nH=%d, nMg=%d] does not have a " 
                            "dissociation table"
                            % (name, cid, nH, nMg))

        diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
        diss_table.SetCharge(nH, charge, nMg)
        dG0_prime = diss_table.Transform(pH, I, pMg, T)
        csv_out.writerow([cid, kegg.cid2name(cid), "%.1f" % dG0_prime, pH, I, pMg, T, True, source, None])
        
if __name__ == "__main__":
    main()