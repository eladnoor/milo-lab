import yeast_prepare_db, pylab
from yeast_stoichiometric_lp import YeastStoichiometricLP

if (__name__ == "__main__"):
    comm = yeast_prepare_db.connect_db()
    cursor = comm.cursor()
    (species, reactions) = yeast_prepare_db.get_model(cursor)
    r_biomass = [r[1] for r in reactions].index('biomass production')
    
    max_flux = 100
    slip = YeastStoichiometricLP("Yeast", self.LOG_FILE)
    slip.add_flux_constraint_list(reactions, add_indicator_variables=True)
        
    
    slip.add_stoichiometric_constraints(f, S, [s[0] for s in species], \
                                        [(r[0], '=>') for r in reactions], \
                                        source=None, target=None)
    slip.set_objective()
    slip.export("../res/pathologic/%s/lp_minflux.txt" % experiment_name)
    print slip.solve()