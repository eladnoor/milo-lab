# This file reads the data stored in the gibbs.sqlite database and 

import os, sys
orig_dir = os.getcwd()
pygibbs_path, _ = os.path.split(orig_dir)
src_path, _ = os.path.split(pygibbs_path)
os.chdir(src_path)
print src_path
sys.path.append(src_path)

import numpy as np
from pygibbs.unified_group_contribution import UnifiedGroupContribution
from toolbox.database import SqliteDatabase

db = SqliteDatabase('../res/gibbs.sqlite', 'r')
ugc = UnifiedGroupContribution(db)
ugc.LoadGroups(FromDatabase=True)
ugc.LoadObservations(FromDatabase=True)
ugc.LoadGroupVectors(FromDatabase=True)
ugc.LoadData(FromDatabase=True)
ugc.init()
result_dict = ugc._GetContributionData(ugc.S.copy(), ugc.cids,
                                       ugc.b.copy(), ugc.anchored)

g_pgc = result_dict['group_contributions']
P_L_pgc = result_dict['pgc_conservations']

os.chdir(orig_dir)
print os.getcwd()
np.save('g_pgc.gz', g_pgc)
np.save('P_L_g_pgc.gz', P_L_pgc)
