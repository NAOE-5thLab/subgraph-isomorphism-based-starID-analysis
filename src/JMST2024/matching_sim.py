import os
import time
import itertools
import multiprocessing

import numpy as np
from scipy.stats import special_ortho_group

from subgraph import (
    YaleStarCatalog,
    StarDatabase,
    PairDatabase,
    SubgraphIsomorphismBasedMatching,
)
from utils import (
    limit_uniform_spherical_vector,
    rodrigues_rotation_matrix,
    equatorial2vec,
)


### probrem setting
# os
log_dir = "./log/matching/"
seed = 100
multi = True
# param
U = 1024
# sampling
N_loop = int(1e1) + 1
# N_loop = int(1e4) + 1
# pattern
theta_FOV_list = [5, 10, 20, 40, 80]
M_lim_list = [3.5, 4.5, 5.5]
beta_list = [0.0, 0.2, 0.4, 0.6, 0.8]
M_lim_max = max(M_lim_list)


def matching_sim(theta_FOV, M_lim, beta):
    path = f"{log_dir}/matching_FOV{theta_FOV}_M{M_lim}_beta{beta}.dat"
    cash_dir = f"{log_dir}/cash/FOV{theta_FOV}_M{M_lim}_beta{beta}/"
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
    catalog = YaleStarCatalog(log_dir=cash_dir)
    df_D_C = catalog.get_df()
    s_vec = equatorial2vec(catalog.get_RA(), catalog.get_DE())
    df_D_C.loc[:, ("s_X", "s_Y", "s_Z")] = s_vec
    # StarDB
    D_DB = StarDatabase(df_D_C)
    D_DB.filtering_by_visual_magnitude(M_lim_max)
    D_DB.filtering_by_multiple_stars(theta_min)
    D_DB_HR = D_DB.get_HR()
    # PairDB
    P_DB = PairDatabase(D_DB.get_df(), log_dir=cash_dir)
    P_DB.create_catalog(theta_max)
    # StarDB (obs)
    D_DB_OBS = StarDatabase(df_D_C)
    D_DB_OBS.filtering_by_visual_magnitude(M_lim_max)
    D_DB_OBS.filtering_by_multiple_stars(theta_min)
    D_DB_OBS.filtering_by_visual_magnitude(M_lim)
    # matching
    matching = SubgraphIsomorphismBasedMatching(D_DB, P_DB, epsilon)

    ### Simulation ###
    with open(path, "w") as f:
        header = ["index", "N_obs", "N_match", "N_candi", "unique", "included"]
        for data in header[:-1]:
            f.write(f"{data},")
        f.write(f"{header[-1]}\n")
    #
    log_list = []
    for i in range(N_loop):
        print(f"\r{os.path.basename(path)[:-3]} : {i}/{N_loop}", end="")
        ### Observed stars ###
        s_vec_hat, obs_HR = get_random_stars(
            D_DB_OBS, np_random, theta_FOV, beta, epsilon
        )
        ### Matching ###
        candi_IDs, obs_IDs, info = matching.match_stars(
            s_vec_hat, np_random, with_check=False
        )
        #
        N_obs = len(s_vec_hat)
        N_match = len(info["N_candi"])
        N_candi = len(candi_IDs)
        #
        unique = N_candi == 1
        included = []
        for candi_ID in candi_IDs:
            included.append(set(obs_HR[obs_IDs]) == set(D_DB_HR[candi_ID]))
        included = any(included)
        ### Save result ###
        log_list.append([i, N_obs, N_match, N_candi, unique, included])
        # save
        freq = 10
        if (i % freq == 0) and (i != 0):
            with open(path, "a") as f:
                for log in log_list[i - freq + 1 : i + 1]:
                    for data in log[:-1]:
                        f.write(f"{data},")
                    f.write(f"{log[-1]}\n")
    return log_list


def get_random_stars(
    D_DB: StarDatabase,
    np_random: np.random.Generator,
    theta_FOV,
    beta,
    epsilon,
):
    s_vec = D_DB.get_df()[["s_X", "s_Y", "s_Z"]].to_numpy()
    N_obs = len(s_vec)
    ### Observed stars ###
    # random rotation
    R = special_ortho_group.rvs(dim=3, size=1, random_state=np_random)
    R_s_vec = np.dot(R, s_vec.T).T
    # FOV condition
    cond_1 = R_s_vec[:, 2] > 0.0
    temp = np.tan(theta_FOV / 2.0) * R_s_vec[:, 2]
    cond_2 = np.abs(R_s_vec[:, 0]) < temp
    cond_3 = np.abs(R_s_vec[:, 1]) < temp
    cond = (cond_1 * cond_2 * cond_3).astype(np.int32)
    # cover rate
    cond *= np_random.binomial(size=N_obs, n=1, p=1 - beta)
    # select stars
    cond = cond == 1
    obs_HR = D_DB.get_HR()[cond]
    obs_RA = D_DB.get_RA()[cond]
    obs_DE = D_DB.get_DE()[cond]
    ### Add noise ###
    obs_noise = limit_uniform_spherical_vector(epsilon / 2, len(obs_HR), np_random)
    axis_vec = np.concatenate(
        [
            np.cos(obs_RA + 0.5 * np.pi)[:, np.newaxis],
            np.sin(obs_RA + 0.5 * np.pi)[:, np.newaxis],
            np.zeros_like(obs_HR)[:, np.newaxis],
        ],
        axis=1,
    )
    R_s_vec = rodrigues_rotation_matrix(0.5 * np.pi - obs_DE, axis_vec)
    if len(obs_HR) == 1:
        s_vec_hat = rotate_(R_s_vec, obs_noise[np.newaxis, :])
    else:
        s_vec_hat = rotate_(R_s_vec, obs_noise)
    return s_vec_hat, obs_HR


def rotate_(R, vec):
    ret = np.empty((vec.shape[0], vec.shape[1]))
    for i in range(R.shape[0]):
        ret[i, :] = np.dot(R[i], vec[i].T).T
    return ret


if __name__ == "__main__":
    patterns = list(itertools.product(theta_FOV_list, M_lim_list, beta_list))
    i = 0
    if multi:
        processes = []
        while True:
            if len(processes) >= multiprocessing.cpu_count() - 1:
                for process in processes:
                    if not process.is_alive():
                        process.join()
                        processes.remove(process)
            else:
                theta_FOV, M_lim, beta = patterns[i]
                i += 1
                #
                process = multiprocessing.Process(
                    target=matching_sim, args=(theta_FOV, M_lim, beta)
                )
                processes.append(process)
                process.start()
                #
                if len(patterns) <= i:
                    break
            time.sleep(1)
    else:
        for theta_FOV, M_lim, beta_list in patterns:
            log_list = matching_sim(theta_FOV, M_lim, beta_list)
    print("Task complete")
