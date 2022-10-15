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
    N_candi = N_pairset
    # 
    N_candi_setid, candi_setid, time = subgraph_f.subgraph_isomorphism_starid_detector.matching_set_for_analysis(
        n_obs=n_obs, s_hat_set=s_hat_set, epsilon=epsilon,
        n_candi=N_candi, n_set=N_set, ra_set=RA_set, de_set=DE_set,
        n_pairset=N_pairset, thetaset=thetaset, pairidset=pairidset)
    #
    candi_setid_each_list = []
    for i in range(n_obs-1):
        candi_setid_each_list.append(candi_setid[i, :N_candi_setid[i], :i+2])
    time_list = []
    for i in range(n_obs):
        time_list.append(time[i])
    return candi_setid_each_list, time_list


# def matching_set(n_obs, s_hat_list, epsilon, db):
#     s_hat_set = np.array(s_hat_list)
#     N_set = len(db.get_I())
#     RA_set = db.get_RA()
#     DE_set = db.get_DE()
#     N_pairset = len(db.get_I_pair_FOV())
#     thetaset = db.get_Theta_pair_FOV()
#     pairidset = db.get_I_pair_FOV()
#     #
#     N_candi_setid, candi_setid = subgraph_f.subgraph_isomorphism_starid_detector.matching_set(
#         n_obs=n_obs, s_hat_set=s_hat_set, epsilon=epsilon,
#         N_set=N_set, RA_set=RA_set, DE_set=DE_set,
#         N_pairset=N_pairset, thetaset=thetaset, pairidset=pairidset)
#     return candi_setid[:N_candi_setid, :]
