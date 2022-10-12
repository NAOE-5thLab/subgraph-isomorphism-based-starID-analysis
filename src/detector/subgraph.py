import numpy as np
import subgraph_f


def matching_set_for_analysis(n_obs, s_hat_list, epsilon, db):
    s_hat_set = np.array(s_hat_list)
    N_set = len(db.get_I())
    RA_set = db.get_RA()
    DE_set = db.get_DE()
    N_pairset = len(db.get_I_pair_FOV())
    thetaset = db.get_Theta_pair_FOV()
    pairidset = db.get_I_pair_FOV()
    #
    N_candi_setid, candi_setid = subgraph_f.matching_set_for_analysis(
        N_candi_setid, candi_setid,
        n_obs, s_hat_set, epsilon,
        N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
    #
    candi_setid_each_list = []
    for n in range(n_obs-1):
        candi_setid_each_list.append(candi_setid[n, :N_candi_setid[n], :n+2])
    return candi_setid_each_list


def matching_set(n_obs, s_hat_list, epsilon, db):
    s_hat_set = np.array(s_hat_list)
    N_set = len(db.get_I())
    RA_set = db.get_RA()
    DE_set = db.get_DE()
    N_pairset = len(db.get_I_pair_FOV())
    thetaset = db.get_Theta_pair_FOV()
    pairidset = db.get_I_pair_FOV()
    #
    N_candi_setid, candi_setid = subgraph_f.matching_set(
        N_candi_setid, candi_setid,
        n_obs, s_hat_set, epsilon,
        N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
    return candi_setid[:N_candi_setid, :]
