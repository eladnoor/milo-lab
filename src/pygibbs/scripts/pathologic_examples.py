import logging
from pygibbs.kegg_reaction import Reaction
from pygibbs.pathologic import Pathologic
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
import sys
from pygibbs.thermodynamic_estimators import LoadAllEstimators


def get_thermo(set_bounds=True):
    estimators = LoadAllEstimators()
    thermo = estimators['merged']
    thermo.c_range = (1e-10, 1e-2)

    if set_bounds:
        thermo.bounds[1]    = (1,    1)    # water
        thermo.bounds[11]   = (1e-5, 1e-5) # CO2
        thermo.bounds[288]  = (9e-5, 9e-5) # carbonate
        thermo.bounds[2]    = (5e-3, 5e-3) # ATP (in order to keep ATP -> ADP at -55 kJ/mol)
        thermo.bounds[8]    = (5e-4, 5e-4) # ADP (in order to keep ATP -> ADP at -55 kJ/mol)
        thermo.bounds[9]    = (5e-3, 5e-3) # Pi  (in order to keep ATP -> ADP at -55 kJ/mol)
        thermo.bounds[20]   = (1e-4, 1e-4) # AMP (in order to keep ATP -> AMP at -110 kJ/mol)
        thermo.bounds[13]   = (3e-9, 3e-9) # PPi (in order to keep ATP -> AMP at -110 kJ/mol)
        #thermo.bounds[20]   = (5e-3, 5e-3), # AMP (in order to keep ATP -> AMP at -55 kJ/mol)
        #thermo.bounds[13]   = (5e-2, 5e-2), # PPi (in order to keep ATP -> AMP at -55 kJ/mol)
        thermo.bounds[3]    = (1e-4, 1e-3) # NAD (ox)
        thermo.bounds[4]    = (1e-4, 1e-3) # NAD (red)
        thermo.bounds[6]    = (1e-4, 1e-3) # NADP (ox)
        thermo.bounds[5]    = (1e-4, 1e-3) # NADP (red)
        thermo.bounds[399]  = (1e-4, 1e-3) # ubiquinone (ox)
        thermo.bounds[390]  = (1e-4, 1e-3) # ubiquinol (red)
        thermo.bounds[139]  = (1e-4, 1e-3) # ferredoxin (ox)
        thermo.bounds[138]  = (1e-4, 1e-3) # ferredoxin (red)
                      
        thermo.bounds[828]  = (1e-4, 1e-3) # menaquinone (ox)
        thermo.bounds[5819] = (1e-4, 1e-3) # menaquinone (red)
        thermo.bounds[343]  = (1e-4, 1e-3) # thioredoxin (ox)
        thermo.bounds[342]  = (1e-4, 1e-3) # thioredoxin (red)
        thermo.bounds[876]  = (1e-4, 1e-3) # coenzyme F420 (ox)
        thermo.bounds[1080] = (1e-4, 1e-3) # coenzyme F420 (red)
    
    return thermo

def example_lower_glycolysis():
    
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    thermo=get_thermo(False),
                    max_solutions=None,
                    max_reactions=8,
                    maximal_dG=0.0,
                    thermodynamic_method='global',
                    update_file=None)
    pl.add_cofactor_reaction(Reaction.FromFormula("C00001 <=> null", name='Free H2O'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00009 <=> null", name='Free Pi'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00013 <=> null", name='Free PPi'))

    # all electron transfer reactions
    pl.add_cofactor_reaction(Reaction.FromFormula("C00003 <=> C00004", name='NAD redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00006 <=> C00005", name='NADP redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00016 <=> C01352", name='FAD redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00138 <=> C00139", name='ferredoxin redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00030 <=> C00028", name='acceptor/donor redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00125 <=> C00126", name='ferricytochrome c redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00996 <=> C00999", name='ferricytochrome b5 redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C01070 <=> C01071", name='ferricytochrome c-553 redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C05906 <=> C01617", name='leucocyanidin redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00343 <=> C00342", name='thioredoxin disulfide redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C03648 <=> C00974", name='cis-3,4-Leucopelargonidin redox'))

    # all Adenosine phosphorylations
    pl.add_cofactor_reaction(Reaction.FromFormula("C00002 <=> C00008", name='ATP to ADP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00002 <=> C00020", name='ATP to AMP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00008 <=> C00020", name='ATP to AMP'))

    pl.add_cofactor_reaction(Reaction.FromFormula("C00131 <=> C00206", name='dATP to dADP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00131 <=> C00360", name='dATP to dAMP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00206 <=> C00360", name='dATP to dAMP'))

    pl.add_cofactor_reaction(Reaction.FromFormula("C00081 <=> C00104", name='ITP to IDP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00081 <=> C00130", name='ITP to IMP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00104 <=> C00130", name='ITP to IMP'))

    pl.add_cofactor_reaction(Reaction.FromFormula("C00044 <=> C00035", name='GTP to GDP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00044 <=> C00144", name='GTP to GMP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00035 <=> C00144", name='GTP to GMP'))

    pl.add_cofactor_reaction(Reaction.FromFormula("C00063 <=> C00112", name='CTP to CDP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00063 <=> C00055", name='CTP to CMP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00112 <=> C00055", name='CTP to CMP'))

    r = Reaction.FromFormula("C00003 + C00118 + C00001 => C00022 + C00004 + C00009")
    #r.Balance()
    pl.find_path("GAP => PYR", r)

def example_oxidative():
    
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    thermo=get_thermo(),
                    max_solutions=None,
                    max_reactions=20,
                    maximal_dG=0.0,
                    thermodynamic_method='global',
                    update_file=None)
    pl.add_cofactor_reaction(Reaction.FromFormula("C00001 <=> null", name='Free H2O'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00002 <=> null", name='Free ATP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00003 <=> null", name='Free NAD+'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00004 <=> null", name='Free NADH'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00005 <=> null", name='Free NADPH'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00006 <=> null", name='Free NADP+'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00008 <=> null", name='Free ADP'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00009 <=> null", name='Free Pi'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00013 <=> null", name='Free PPi'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00020 <=> null", name='Free AMP'))
    
    if False:
        pl.add_cofactor_reaction(Reaction.FromFormula("C00009 <=> null", name='Free Pi'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00013 <=> null", name='Free PPi'))
    
        # all electron transfer reactions
        pl.add_cofactor_reaction(Reaction.FromFormula("C00003 <=> C00004", name='NAD redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00006 <=> C00005", name='NADP redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00016 <=> C01352", name='FAD redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00138 <=> C00139", name='ferredoxin redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00030 <=> C00028", name='acceptor/donor redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00125 <=> C00126", name='ferricytochrome c redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00996 <=> C00999", name='ferricytochrome b5 redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C01070 <=> C01071", name='ferricytochrome c-553 redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C05906 <=> C01617", name='leucocyanidin redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00343 <=> C00342", name='thioredoxin disulfide redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C03648 <=> C00974", name='cis-3,4-Leucopelargonidin redox'))
    
        # all Adenosine phosphorylations
        pl.add_cofactor_reaction(Reaction.FromFormula("C00002 <=> C00008", name='ATP to ADP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00002 <=> C00020", name='ATP to AMP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00008 <=> C00020", name='ATP to AMP'))
    
        pl.add_cofactor_reaction(Reaction.FromFormula("C00131 <=> C00206", name='dATP to dADP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00131 <=> C00360", name='dATP to dAMP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00206 <=> C00360", name='dATP to dAMP'))
    
        pl.add_cofactor_reaction(Reaction.FromFormula("C00081 <=> C00104", name='ITP to IDP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00081 <=> C00130", name='ITP to IMP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00104 <=> C00130", name='ITP to IMP'))
    
        pl.add_cofactor_reaction(Reaction.FromFormula("C00044 <=> C00035", name='GTP to GDP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00044 <=> C00144", name='GTP to GMP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00035 <=> C00144", name='GTP to GMP'))
    
        pl.add_cofactor_reaction(Reaction.FromFormula("C00063 <=> C00112", name='CTP to CDP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00063 <=> C00055", name='CTP to CMP'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00112 <=> C00055", name='CTP to CMP'))

    #r = Reaction.FromFormula("C00003 + C00118 + C00001 => C00022 + C00004 + C00009")
    #r = Reaction.FromFormula("2 C00022 => C00092")
    r = Reaction.FromFormula("C00092 => 3 C00011")
    #r.Balance()
    pl.find_path("oxidative", r)

def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    #example_lower_glycolysis()
    example_oxidative()

if __name__ == '__main__':
    main()

# Handy reference
#
# name         =  CID
# ---------------------
# atp          = C00002
# adp          = C00008
# pi           = C00009
# co2          = C00011
# nad          = C00003
# nadh         = C00004
# glucose      = C00031
# g6p          = C00092
# fbp          = C00354
# bpg          = C00236
# g3p          = C00118
# threepg      = C00197
# pep          = C00074
# pyruvate     = C00022
# succinyl_coa = C00091
# acetyl_coa   = C00024
# lactate      = C00186
# acetate      = C00033
# carbon       = C19202
