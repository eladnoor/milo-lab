from groups import GroupContribution

G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="dG0_train")
G.load_groups("../data/groups_species.csv")
G.train("../data/dG0.csv", use_dG0_format=True)
G.analyze_training_set()
G.load_hatzimanikatis_rid_data("../data/hatzimanikatis_rid.csv")
G.load_hatzimanikatis_cid_data("../data/hatzimanikatis_cid.csv")
