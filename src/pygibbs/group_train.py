#!/usr/bin/python

from pygibbs.groups import GroupContribution
from toolbox.sql import write_csv_table

G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_train")
G.load_groups("../data/thermodynamics/groups_species.csv")
G.train("../data/thermodynamics/dG0.csv", use_dG0_format=True)
G.analyze_training_set()
G.load_cid2pmap(recalculate=True)
#write_csv_table(G.comm.cursor(), "../res/dG0_train/estimations.csv", "gc_cid2prm")
#write_csv_table(G.comm.cursor(), "../res/dG0_train/errors.csv", "gc_cid2error")
