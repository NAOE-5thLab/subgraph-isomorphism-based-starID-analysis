import os
import argparse
import tqdm

import numpy as np
import pandas as pd

import log
import catalog
import db
import rand
from config import Param
from detector.subgraph import matching_set_for_analysis
from utils.seeds import gen_seeds
from utils.equatorial import equatorial2vec
from utils.interangle import inter_star_angle_vec


def calc_stats(candi_setids, obs_setid):
    obs_set = set(obs_setid)
    candi_set_list = [set(candi_setid) for candi_setid in candi_setids]
    # matching number
    matching_num = len(candi_setids)
    multiple = matching_num > 1
    unique = matching_num == 1
    noexist = matching_num == 0
    # correct pair is included in candidate pairs
    included = obs_set in candi_set_list
    #
    if len(candi_set_list) == 0:
        candi_intersect_set = set()
    else:
        candi_intersect_set = set.intersection(*candi_set_list)
    #
    if candi_intersect_set == set():
        determined = False
        correct = False
    else:
        determined = True
        correct = candi_intersect_set <= obs_set
    return matching_num, multiple, unique, noexist, included, determined, correct


def main(args):
    ### Prepair PARAMETER ###
    conf = Param()
    Vmax = conf.Vmax_list[args.i_Vmax]
    theta_FOV = conf.theta_FOV_list[args.i_theta_FOV]
    theta_img = conf.theta_img_list[args.i_theta_img]
    k = conf.k_list[args.i_k]
    kappa = rand.estimate_kappa_von_mises_fisher_3d(theta_img, conf.alpha)

    i_parallel = args.i_parallel

    ### Prepair DATABASE ###
    # origne catalog
    yale_catalog = catalog.YaleStarCatalog()
    # star db
    star_db = db.StarDB(
        yale_catalog.get_HR(), yale_catalog.get_RA(),
        yale_catalog.get_DE(), yale_catalog.get_Vmag(), Vmax=Vmax)
    I = star_db.get_I()
    RA = star_db.get_RA()
    DE = star_db.get_DE()
    # pair star db
    pairstar_db = db.PairStarDB(I, RA, DE, theta_FOV=theta_FOV)

    ### Prepair LOGGER ###
    header = [
        'seed', 'observable', 'time', 'matching_num',
        'multiple', 'unique', 'noexist',
        'included', 'determined', 'correct'
    ]
    logger_list = []
    for i, n_obs in enumerate(conf.n_obs_list):
        logger = log.Logger(conf.log_temp_dir, header)
        logger.reset()
        logger_list.append(logger)

    ### MONTE CARLO SIMULATION ###
    calc_num = int(len(conf.calc_seeds)/conf.parallel_num)
    if i_parallel + 1 == conf.parallel_num:
        calc_seeds = conf.calc_seeds[i_parallel*calc_num:]
    else:
        calc_seeds = conf.calc_seeds[
            i_parallel*calc_num:(i_parallel+1)*calc_num]
    #
    for sim_seed in tqdm.tqdm(calc_seeds):
        ### CREATE SEED ###
        sim_seeds = gen_seeds(sim_seed, 3)
        ### RANDOM SAMPLING ###
        # select center randomly
        center_vec = rand.uniform_spherical_vector(seed=sim_seeds[0])
        # collect stars in circle (FOV)
        in_circle = inter_star_angle_vec(
            center_vec, equatorial2vec(RA, DE)) < theta_FOV
        ID_in_circle_list = list(I[in_circle])
        # select stars randomly
        obs_setid = rand.random_select(
            conf.n_obs_list[-1], ID_in_circle_list, seed=sim_seeds[1])
        n_obs_able = len(obs_setid)
        s_list = list(equatorial2vec(RA[obs_setid], DE[obs_setid]))
        # add noise on star position
        s_hat_list = []
        theta_error_list = []
        noise_seeds = gen_seeds(sim_seeds[2], len(s_list))
        for i, s in enumerate(s_list):
            s_hat, theta_error = rand.von_mises_fisher_3d_sampling(
                s, kappa, seed=noise_seeds[i], theta_flag=True)
            s_hat_list.append(s_hat)
            theta_error_list.append(theta_error)
        ### SEARCH GRAPH ###
        if n_obs_able >= 2:
            # matching
            candi_setid_each_list, time_list = matching_set_for_analysis(
                n_obs_able, s_hat_list, k*theta_img, pairstar_db)
        #
        for i, n_obs in enumerate(conf.n_obs_list):
            if n_obs <= n_obs_able:
                # stats
                matching_num, multiple, unique, noexist, included, determined, correct = calc_stats(
                    candi_setid_each_list[i], obs_setid[:i+2])
                if i == 0:
                    time = time_list[1] - time_list[0]
                else:
                    time = time_list[i+1] - time_list[1]
                # log
                logger_list[i].append(
                    [sim_seed, 1, time, matching_num, int(multiple), int(unique),
                     int(noexist), int(included), int(determined), int(correct)])
            else:
                logger_list[i].append(
                    [sim_seed, 0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0])
    # save
    for i, n_obs in enumerate(conf.n_obs_list):
        fname = f'stats_{i_parallel}_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}'
        logger_list[i].save(fname)

    ### Marge Results ###
    finish = True
    for i in range(conf.parallel_num):
        fname = f'stats_{i}_{conf.n_obs_list[-1]}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}'
        if not os.path.isfile(conf.log_temp_dir + '/' + f'{fname}.csv'):
            finish = False
            break
    if finish:
        for n_obs in conf.n_obs_list:
            data_list = []
            for i in range(conf.parallel_num):
                fname = f'stats_{i}_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}'
                data = pd.read_csv(
                    conf.log_temp_dir + '/' + f'{fname}.csv', header=0, index_col=0).to_numpy()
                data_list.append(data)
                #
                os.remove(conf.log_temp_dir + '/' + f'{fname}.csv')
            marged_data = np.concatenate(data_list, axis=0)
            # save
            df = pd.DataFrame(marged_data, columns=header)
            marged_fname = f'stats_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}'
            df.to_csv(conf.log_dir + marged_fname + '.csv')


if __name__ == '__main__':
    # argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--i_Vmax",
        type=int, default=4,
        help="ganbare")
    parser.add_argument(
        "--i_theta_FOV",
        type=int, default=2,
        help="ganbare")
    parser.add_argument(
        "--i_theta_img",
        type=int, default=4,
        help="ganbare")
    parser.add_argument(
        "--i_k",
        type=int, default=3,
        help="ganbare")
    parser.add_argument(
        "--i_parallel",
        type=int, default=0,
        help="ganbare")
    args = parser.parse_args()
    main(args)
    print('Task complete')
