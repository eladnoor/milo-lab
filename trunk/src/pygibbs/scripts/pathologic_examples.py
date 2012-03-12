import logging
from pygibbs.kegg_reaction import Reaction
from pygibbs.pathologic import Pathologic
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
import sys
from pygibbs.thermodynamic_estimators import LoadAllEstimators

def add_cofactor_reactions(pl, NAD_only=False):
    pl.add_cofactor_reaction(Reaction.FromFormula("C00001 <=> null", name='Free H2O'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00009 <=> null", name='Free Pi'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00013 <=> null", name='Free PPi'))

    # all electron transfer reactions
    pl.add_cofactor_reaction(Reaction.FromFormula("C00003 <=> C00004", name='NAD redox'))
    pl.add_cofactor_reaction(Reaction.FromFormula("C00006 <=> C00005", name='NADP redox'))
    if not NAD_only:
        pl.add_cofactor_reaction(Reaction.FromFormula("C00016 <=> C01352", name='FAD redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00138 <=> C00139", name='ferredoxin redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00030 <=> C00028", name='acceptor/donor redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00125 <=> C00126", name='ferricytochrome c redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00996 <=> C00999", name='ferricytochrome b5 redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C01070 <=> C01071", name='ferricytochrome c-553 redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C05906 <=> C01617", name='leucocyanidin redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C00343 <=> C00342", name='thioredoxin disulfide redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C03648 <=> C00974", name='cis-3,4-Leucopelargonidin redox'))
        pl.add_cofactor_reaction(Reaction.FromFormula("C05684 <=> C01528", name='selenide redox'))
    else:
        pl.ban_compound(16)    # FAD(ox)
        pl.ban_compound(1352)  # FAD(red)
        pl.ban_compound(138)   # ferredoxin (ox)
        pl.ban_compound(139)   # ferredoxin (red)
        pl.ban_compound(30)    # acceptor
        pl.ban_compound(28)    # donor
        pl.ban_compound(125)   # ferricytochrome c (ox)
        pl.ban_compound(126)   # ferricytochrome c (red)
        pl.ban_compound(996)   # ferricytochrome b5 (ox)
        pl.ban_compound(999)   # ferricytochrome b5 (red)
        pl.ban_compound(1070)  # ferricytochrome c-553 (ox)
        pl.ban_compound(1071)  # ferricytochrome c-553 (red)
        pl.ban_compound(5906)  # leucocyanidin (ox)
        pl.ban_compound(1617)  # leucocyanidin (red)
        pl.ban_compound(343)   # thioredoxin disulfide (ox)
        pl.ban_compound(342)   # thioredoxin disulfide (red)
        pl.ban_compound(3648)  # cis-3,4-Leucopelargonidin (ox)
        pl.ban_compound(974)   # cis-3,4-Leucopelargonidin (red)
        pl.ban_compound(5684)  # selenide (ox)
        pl.ban_compound(1528)  # selenide (red)

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

def example_lower_glycolysis(thermo):
    
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    thermo=thermo,
                    max_solutions=None,
                    max_reactions=8,
                    maximal_dG=0.0,
                    thermodynamic_method='global',
                    update_file=None)
    add_cofactor_reactions(pl)
    #r = Reaction.FromFormula("C00003 + C00118 + C00001 => C00022 + C00004 + C00009")
    r = Reaction.FromFormula("C00118 => C00022")
    #r.Balance()
    pl.find_path("GAP => PYR", r)

def example_oxidative(thermo):
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    thermo=thermo,
                    max_solutions=None,
                    max_reactions=20,
                    maximal_dG=0,
                    thermodynamic_method='global',
                    update_file=None)
    add_cofactor_reactions(pl, NAD_only=False)
    r = Reaction.FromFormula("C00022 => 3 C00011")
    #r.Balance()
    pl.find_path("oxidative", r)

def example_reductive(thermo):
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    thermo=thermo,
                    max_solutions=None,
                    max_reactions=15,
                    maximal_dG=0.0,
                    thermodynamic_method='global',
                    update_file=None)
    add_cofactor_reactions(pl)
    r = Reaction.FromFormula("3 C00011 => C00022")
    #r.Balance()
    pl.find_path("reductive", r)

def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    estimators = LoadAllEstimators()
    thermo = estimators['UGC']
    thermo.SetConditions(I=0.1)
    #example_lower_glycolysis(thermo)
    example_oxidative(thermo)

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
