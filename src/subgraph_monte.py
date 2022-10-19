import argparse
import tqdm

import numpy as np

import log
import catalog
import db
import rand
from detector.subgraph import matching_set_for_analysis
from utils.seeds import gen_seeds
from utils.equatorial import equatorial2vec
from utils.interangle import inter_star_angle_vec


sampling_type = 0


class Param:
    ### Hyperparameter ###
    # DB
    Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
    theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
    # simulation
    theta_img_list = [np.pi/180*10**i for i in [
        -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
    # subgraph matching
    k_list = [2.0**i for i in [-1.0, -0.5, 0.0, 0.5, 1.0]]
    n_obs_max = 6

    ### basic parameter ###
    # system
    log_dir = './log/subgraph_monte/'
    seed = 10
    # simulation
    sample_N = 10000
    alpha = 0.99
    alpha_sampling = 0.9999

    def __init__(self, args):
        ### Hyperparameter ###
        self.Vmax = self.Vmax_list[args.Vmax]
        self.theta_FOV = self.theta_FOV_list[args.theta_FOV]
        self.theta_img = self.theta_img_list[args.theta_img]
        self.k = self.k_list[args.k]
        ### basic parameter ###
        self.kappa = rand.estimate_kappa_von_mises_fisher_3d(
            self.theta_img, self.alpha)
        self.theta_sampling = rand.estimate_thata_von_mises_fisher_3d(
            self.kappa, self.alpha_sampling)
        self.kappa_sampling = rand.estimate_kappa_von_mises_fisher_3d(
            self.theta_img, self.alpha_sampling)


def add_observation_noise(s_list, kappa, theta_sampling, kappa_sampling, seed=100):
    seeds_noise = gen_seeds(seed, len(s_list))
    #
    s_hat_list = []
    theta_error_list = []
    for i, s in enumerate(s_list):
        if sampling_type == 0:
            s_hat, theta_error = rand.von_mises_fisher_3d_sampling(
                s, kappa, seed=seeds_noise[i], theta_flag=True)
        elif sampling_type == 1:
            s_hat, theta_error = rand.uniform_limited_spherical_vector(
                s, theta_sampling, seed=seeds_noise[i], theta_flag=True)
        elif sampling_type == 2:
            s_hat, theta_error = rand.von_mises_fisher_3d_sampling(
                s, kappa_sampling, seed=seeds_noise[i], theta_flag=True)
        s_hat_list.append(s_hat)
        theta_error_list.append(theta_error)
    return s_hat_list, theta_error_list


def calc_weight(theta, kappa, theta_sampling, kappa_sampling):
    if sampling_type == 0:
        return np.ones_like(theta)
    elif sampling_type == 1:
        numerator = kappa * \
            np.exp(kappa*(np.cos(theta)-1)) * (1 - np.cos(theta_sampling))
        denominator = 1 - np.exp(-2*kappa)
        return numerator/denominator
    elif sampling_type == 2:
        numerator = kappa * (1 - np.exp(-2*kappa_sampling))
        denominator = kappa_sampling * (1 - np.exp(-2*kappa))
        return (numerator/denominator) * np.exp((kappa - kappa_sampling)*(np.cos(theta) - 1))


def main(args):
    ### Prepair PARAMETER ###
    p = Param(args)

    ### Prepair DATABASE ###
    # origne catalog
    yale_catalog = catalog.YaleStarCatalog()
    # star db
    star_db = db.StarDB(
        yale_catalog.get_HR(), yale_catalog.get_RA(),
        yale_catalog.get_DE(), yale_catalog.get_Vmag(), Vmax=p.Vmax)
    I = star_db.get_I()
    RA = star_db.get_RA()
    DE = star_db.get_DE()
    # pair star db
    pairstar_db = db.PairStarDB(
        I, RA, DE, theta_FOV=p.theta_FOV)

    ### Prepair LOGGER ###
    header = [
        'observable', 'time', 'matching_num', 'multiple', 'unique', 'noexist', 'included', 'weight']
    logger_list = []
    for i in range(p.n_obs_max-1):
        logger = log.Logger(p.log_dir, header)
        logger.reset()
        logger_list.append(logger)

    ### MONTE CARLO SIMULATION ###
    # seed
    seeds_monte = gen_seeds(p.seed, p.sample_N)
    for sample_i in tqdm.tqdm(range(p.sample_N)):
        # seed
        seeds_one = gen_seeds(seeds_monte[sample_i], 3)
        ### RANDOM SAMPLING ###
        # select center randomly
        center_vec = rand.uniform_spherical_vector(seed=seeds_one[0])
        # collect stars in circle (FOV)
        in_circle = inter_star_angle_vec(
            center_vec, equatorial2vec(RA, DE)) < p.theta_FOV
        ID_in_circle_list = list(I[in_circle])
        # select stars randomly
        obs_setid = rand.random_select(
            p.n_obs_max, ID_in_circle_list, seed=seeds_one[1])
        n_obs_able = len(obs_setid)
        s_list = list(equatorial2vec(RA[obs_setid], DE[obs_setid]))
        # add noise on star position
        s_hat_list, theta_error_list = add_observation_noise(
            s_list, p.kappa, p.theta_sampling, p.kappa_sampling, seed=seeds_one[2])
        if n_obs_able >= 2:
            ### SEARCH GRAPH ###
            candi_setid_each_list, time_list = matching_set_for_analysis(
                n_obs_able, s_hat_list, p.k*p.theta_img, pairstar_db)
            ### CALC STATS and LOGGING ###
            # weight
            weights = calc_weight(
                np.array(theta_error_list), p.kappa, p.theta_sampling, p.kappa_sampling)
        #
        for i in range(p.n_obs_max-1):
            if i+2 <= n_obs_able:
                # stats
                matching_num, multiple, unique, noexist, included = calc_stats(
                    candi_setid_each_list[i], obs_setid)
                weight = np.prod(weights[:i+2])
                if i == 0:
                    time = time_list[i+1] - time_list[0]
                else:
                    time = time_list[i+1] - time_list[1]
                # log
                logger_list[i].append(
                    [1, time, matching_num, int(multiple), int(unique), int(noexist), int(included), weight])
            else:
                logger_list[i].append([0, -1, -1, -1, -1, -1, -1, -1])
        #
        # if n_obs_able >= 2:
        #     if time_list[-1] - time_list[0] > 60.0:
        #         print('time over')
        #         break
    # save
    for i in range(p.n_obs_max-1):
        logger_list[i].save(
            f'stats_{p.sample_N}_{p.seed}_{i+2}_{p.Vmax}_{p.theta_FOV*180/np.pi}_{p.theta_img*180/np.pi}_{p.k}_{sampling_type}')


def calc_stats(candi_setids, obs_setid):
    # matching number
    matching_num = len(candi_setids)
    multiple = matching_num > 1
    unique = matching_num == 1
    noexist = matching_num == 0
    # correct pair is included in candidate pairs
    included = False
    for candi_setid in candi_setids:
        if set(obs_setid) == set(candi_setid):
            included = True
    return matching_num, multiple, unique, noexist, included


if __name__ == '__main__':
    # argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--Vmax",
        type=int,
        default=4,
        help="ganbare")
    parser.add_argument(
        "--theta_FOV",
        type=int,
        default=2,
        help="ganbare")
    parser.add_argument(
        "--theta_img",
        type=int,
        default=4,
        help="ganbare")
    parser.add_argument(
        "--k",
        type=int,
        default=3,
        help="ganbare")
    args = parser.parse_args()
    #
    main(args)
    print('Task complete')
