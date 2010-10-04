from groups import GroupContribution

G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_kegg")
G.write_gc_tables()
G.init()

pH = [6.0, 6.5, 7.0, 7.5, 8.0]
I = [0.0, 0.1, 0.2]
G.analyze_all_kegg_compounds(pH=pH, I=I)
G.analyze_all_kegg_reactions(pH=pH, I=I)
