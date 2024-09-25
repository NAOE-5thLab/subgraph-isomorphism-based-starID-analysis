import itertools

import numpy as np
from scipy.stats import special_ortho_group
from scipy.stats import binom

from subgraph import YaleStarCatalog, StarDatabase
from utils import *


### probrem setting
# os
log_dir = "./log/obs_stars/"
seed = 100
# param
U = 1024
# sampling
N_loop = int(1e4)
N_batch = int(1e4)
# pattern
theta_FOV_list = [5, 10, 20, 40, 80]
M_lim_list = [3.5, 4.5, 5.5]
beta_list = [0.0, 0.2, 0.4, 0.6, 0.8]
M_lim_max = max(M_lim_list)


def count_sim_ideal(theta_FOV, M_lim):
    print(f"----- count_sim_ideal : theta_FOV={theta_FOV}, M_lim={M_lim} -----")
    theta_FOV = theta_FOV * np.pi / 180
    ### compute params
    seed_seq = np.random.SeedSequence(seed)
    np_random = np.random.Generator(np.random.PCG64(seed_seq))
    theta_res = np.arctan2(2 * np.tan(theta_FOV / 2), U)
    epsilon = 2 * np.sqrt(2) * theta_res
    theta_min = epsilon
    theta_max = 2 * np.arctan(np.sqrt(2) * np.tan(theta_FOV / 2))

    ### Database ###
    # Catalog
    catalog = YaleStarCatalog(log_dir=log_dir)
    df_D_C = catalog.get_df()
    s_vec = equatorial2vec(catalog.get_RA(), catalog.get_DE())
    df_D_C.loc[:, ("s_X", "s_Y", "s_Z")] = s_vec
    # StarDB
    D_DB = StarDatabase(df_D_C)
    D_DB.filtering_by_visual_magnitude(M_lim_max)
    D_DB.filtering_by_multiple_stars(theta_min)
    D_DB.filtering_by_visual_magnitude(M_lim)
    D_DB.get_info()

    ### Simulation ###
    s_vec = D_DB.get_s_vec()
    obs_stars_N = []
    for i in range(N_loop):
        print(f"\r{i}/{N_loop}", end="")
        # random rotation
        R = special_ortho_group.rvs(dim=3, size=N_batch, random_state=np_random)
        R_s_vec = rotate(R, s_vec)
        R_s_vec = np.dot(R, s_vec.T).T
        ez = np.array([0, 0, 1])
        ezR = np.dot(ez, R)
        # FOV condition
        cond_1 = R_s_vec[:, :, 2] > 0.0
        temp = np.tan(theta_FOV / 2.0) * R_s_vec[:, :, 2]
        cond_2 = np.abs(R_s_vec[:, :, 0]) < temp
        cond_3 = np.abs(R_s_vec[:, :, 1]) < temp
        cond = cond_1 * cond_2 * cond_3
        obs_stars_N += np.sum(cond, axis=1).tolist()
    return obs_stars_N


def hist_with_covered(obs_stars_N, beta):
    print(f"===== hist_with_covered : beta={beta} =====")
    # calc hist
    hist_count = np.zeros(10000)
    max_stars_n = 100
    for obs_stars_n in obs_stars_N:
        hist_count[obs_stars_n] += 1
        if max_stars_n < obs_stars_n:
            max_stars_n = obs_stars_n
    #
    pN_on_thetaFOV_Mlim = hist_count / len(obs_stars_N)
    pN_on_thetaFOV_Mlim_beta = np.zeros(max_stars_n)
    for n in range(max_stars_n):
        for N in range(max_stars_n):
            if pN_on_thetaFOV_Mlim[N] > 0:
                pn_N = binom.pmf(n, N, 1 - beta)
                pN_on_thetaFOV_Mlim_beta[n] += pN_on_thetaFOV_Mlim[N] * pn_N
            else:
                pN_on_thetaFOV_Mlim_beta[n] += 0.0
    assert np.abs(pN_on_thetaFOV_Mlim.sum() - 1) < 1e-3
    assert np.abs(pN_on_thetaFOV_Mlim_beta.sum() - 1) < 1e-3
    return pN_on_thetaFOV_Mlim_beta


def rotate(R, vec):
    ret = np.empty((R.shape[0], vec.shape[0], vec.shape[1]))
    for i in range(R.shape[0]):
        ret[i, :, :] = np.dot(R[i], vec.T).T
    return ret


if __name__ == "__main__":
    for theta_FOV, M_lim in list(itertools.product(theta_FOV_list, M_lim_list)):
        obs_stars_N = count_sim_ideal(theta_FOV, M_lim)
        # save
        path = f"{log_dir}/count_FOV{theta_FOV}_M{M_lim}.dat"
        with open(path, "w") as f:
            for obs_stars_n in obs_stars_N:
                f.write(f"{obs_stars_n}\n")
        for beta in beta_list:
            path = f"{log_dir}/hist_FOV{theta_FOV}_M{M_lim}_beta{beta}.dat"
            pN_on_thetaFOV_Mlim_beta = hist_with_covered(obs_stars_N, beta)
            # save
            with open(path, "w") as f:
                for p in pN_on_thetaFOV_Mlim_beta:
                    f.write(f"{p}\n")
    print("task complete")
