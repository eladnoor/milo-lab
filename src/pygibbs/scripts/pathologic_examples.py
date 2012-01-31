import logging
from pygibbs.kegg_reaction import Reaction
from pygibbs.pathologic import Pathologic
from toolbox.database import SqliteDatabase
from toolbox.html_writer import HtmlWriter
import sys

def example_lower_glycolysis():
    pl = Pathologic(db=SqliteDatabase('../res/gibbs.sqlite', 'r'),
                    public_db=SqliteDatabase('../data/public_data.sqlite'),
                    html_writer=HtmlWriter('../res/pathologic.html'),
                    max_solutions=None,
                    max_reactions=6,
                    thermodynamic_method='global',
                    update_file='../data/thermodynamics/database_updates_fermentation.txt')
    #pl.add_reaction(Reaction.FromFormula("C00009 <=> null", name='Free Pi'))
    #pl.add_reaction(Reaction.FromFormula("C00001 <=> null", name='Free H2O'))
    #pl.add_reaction(Reaction.FromFormula("C00003 <=> C00004", name='Free NAD redox'))
    #pl.add_reaction(Reaction.FromFormula("C00006 <=> C00005", name='Free NADP redox'))
    #r = Reaction.FromFormula("C00003 + C00118 + C00001 => C00022 + C00004 + C00009")
    r = Reaction.FromFormula("C00118 => C00022")
    #r.Balance()
    pl.find_path("GA3P => PYR", r)

def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    example_lower_glycolysis()

if __name__ == '__main__':
    main()

# Handy reference
#
# name         =  CID
# ---------------------
# atp          = C00002
# adp          = C00008
# pi           = C00009
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
