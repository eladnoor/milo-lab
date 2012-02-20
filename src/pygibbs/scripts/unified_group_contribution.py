import numpy as np
from pygibbs.group_vector import GroupVector
from pygibbs.groups_data import GroupsData
from pygibbs.group_observation import GroupObervationCollection
from toolbox.linear_regression import LinearRegression
import logging
from toolbox.database import SqliteDatabase

def LoadData(db):
    groups_data = GroupsData.FromDatabase(db, transformed=False)

    cid2nH_nMg = {}
    cid2error = {}
    cid2groupvec = {}
    
    for row in db.DictReader('pgc_groupvector'):
        cid = row['cid']
        gv_str = row['groupvec']
        if gv_str is not None:
            groupvec = GroupVector.FromJSONString(groups_data, gv_str)
        else:
            groupvec = None
        cid2groupvec[cid] = groupvec
        cid2error[cid] = row['err']
        cid2nH_nMg[cid] = (row['nH'], row['nMg'])
    
    obs_collection = GroupObervationCollection.FromDatabase(db=db,
                            table_name='pgc_observations', transformed=False)
    cids, S, b, anchored = obs_collection.GetStoichiometry()

    anchored_cols = np.where(anchored==1)[1]
    # now remove anchored data from S and leave only the data which will be 
    # used for calculating the group contributions
    g, _ = LinearRegression.LeastSquares(S[:, anchored_cols],
                                         b[:, anchored_cols])
    P_C, P_L = LinearRegression.ColumnProjection(S[:, anchored_cols])
    b -= np.dot(np.dot(g, P_C), S)
    S = np.dot(P_L, S)
    
    # set epsilon-small values to absolute 0
    S[np.where(abs(S) < 1e-10)] = 0
    
    # removed zero rows (compounds) from S
    used_cid_indices = set(np.nonzero(np.sum(abs(S), 1))[0])
    for i_cid, cid in enumerate(cids):
        if cid2groupvec[cid] is None:
            used_cid_indices.difference_update([i_cid])
            for i_obs in np.nonzero(S[i_cid, :])[0]:
                S[:, i_obs] = 0

    used_cid_indices = sorted(used_cid_indices)
    S = S[used_cid_indices, :]
    n_groups = len(groups_data.GetGroupNames()) # number of groups
    G = np.zeros((len(used_cid_indices), n_groups))
    for i, i_cid in enumerate(used_cid_indices):
        G[i, :] = cid2groupvec[cids[i_cid]].Flatten()
    
    return G, S, b
    
def Estimate(G, S, b, r):
    """
        Given
            - group matrix: G (m x g)
            - stoichiometric matrix: S (m x n)
            - Gibbs energies: b (1 x n) 
            - a new reaction: r (m x 1)
        Return:
            The estimated Gibbs energy of 'r'
    """
    P_C, P_L = LinearRegression.ColumnProjection(S)
    cont_S, _ = LinearRegression.LeastSquares(S, b)
    cont_GS, ker_GS = LinearRegression.LeastSquares(np.dot(G.T, S), b)
    
    # we must check that r is orthogonal to ker_GS
    
    cont_S = np.dot(cont_S, P_C)
    cont_GS = np.dot(np.dot(cont_GS, G.T), P_L)
    return np.dot(cont_S, r)[0], np.dot(cont_GS, r)[0]

def LOO(G, S, b):
    n = S.shape[1]
    for i in xrange(n):
        print b[0, i]
        no_i = range(0, i) + range(i+1, n)
        print Estimate(G, S[:, no_i], b[:, no_i], S[:, i])
    
if __name__ == "__main__":
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    G, S, b = LoadData(db)
    LOO(G, S, b)
    