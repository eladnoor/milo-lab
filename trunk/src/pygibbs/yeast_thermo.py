import pylab
from pygibbs import yeast_prepare_db
from pygibbs.yeast_stoichiometric_lp import YeastStoichiometricLP

if (__name__ == "__main__"):
    comm = yeast_prepare_db.connect_db()
    cursor = comm.cursor()
    (species, reactions) = yeast_prepare_db.get_model(cursor)
    biomass_rid = yeast_prepare_db.get_name2rid(cursor)['biomass production']
    
    max_flux = 100
    slip = YeastStoichiometricLP("Yeast", milp=True)
    slip.add_flux_constraint_list(reactions)
    slip.add_margin_dGf_constraints(c_mid=1e-4)
    slip.set_pCr_objective()
    slip.export("../res/yeast.lp")
    print slip.solve()