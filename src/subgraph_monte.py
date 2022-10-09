import argparse
import numpy as np

import log
import catalog
import db
import rand
from detector.subgraph import matching_subgraph
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
        -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0]]
    # subgraph matching
    k_list = [2.0**i for i in [-1.0, -0.5, 0.0, 0.5, 1.0]]

    ### basic parameter ###
    # system
    log_dir = './log/'
    seed = 10
    # simulation
    sample_N = 10000
    alpha = 0.99
    alpha_sampling = 0.9999
    # subgraph matching
    n_obs = 2

    def __init__(self, args):
        ### Hyperparameter ###
        self.Vmax = self.Vmax_list[args.Vmax]
        self.theta_FOV = self.theta_FOV_list[args.theta_FOV]
        self.theta_img = self.theta_img_list[args.theta_img]
        self.k = self.k_list[args.k]
        ### basic parameter ###
        self.kappa = rand.estimate_kappa_von_mises_fisher_3d(
            self.theta_img, self.alpha)
        self.kappa_sampling = rand.estimate_kappa_von_mises_fisher_3d(
            self.theta_img, self.alpha_sampling)
        self.theta_sampling = rand.estimate_thata_von_mises_fisher_3d(
            self.kappa, self.alpha_sampling)


def add_observation_noise(s_list, p, seed=100):
    seeds_noise = gen_seeds(seed, p.n_obs)
    s_hat_list = []
    theta_error_list = []
    for i, s in enumerate(s_list):
        if sampling_type == 0:
            s_hat, theta_error = rand.von_mises_fisher_3d_sampling(
                s, p.kappa, seed=seeds_noise[i])
        elif sampling_type == 1:
            s_hat, theta_error = rand.uniform_limited_spherical_vector(
                s, p.theta_sampling, seed=seeds_noise[i], theta_flag=True)
        elif sampling_type == 2:
            s_hat, theta_error = rand.von_mises_fisher_3d_sampling(
                s, p.kappa_sampling, seed=seeds_noise[i])
        s_hat_list.append(s_hat)
        theta_error_list.append(theta_error)


def calc_weight(theta, p):
    if sampling_type == 0:
        return 1
    elif sampling_type == 1:
        numerator = p.kappa * \
            np.exp(p.kappa*(np.cos(theta)-1)) * (1 - np.cos(p.theta_sampling))
        denominator = 1 - np.exp(-2*p.kappa)
        return numerator/denominator
    elif sampling_type == 2:
        numerator = p.kappa * (1 - np.exp(-2*p.kappa_sampling))
        denominator = p.kappa_sampling * (1 - np.exp(-2*p.kappa))
        return (numerator/denominator) * np.exp((p.kappa - p.kappa_sampling)*(np.cos(theta) - 1))


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
        'matching_num', 'multiple', 'unique', 'noexist', 'included', 'weight']
    logger_case_value = log.Logger(p.log_dir, header)
    logger_case_value.reset()

    ### MONTE CARLO SIMULATION ###
    # seed
    seeds_monte = gen_seeds(p.seed, p.sample_N)
    for sample_i in range(p.sample_N):
        seeds_one = gen_seeds(seeds_monte[sample_i], 3)
        ### RANDOM SAMPLING ###
        # select center randomly
        center_vec = rand.uniform_spherical_vector(seed=seeds_one[0])
        # collect stars in circle (FOV)
        s_vec = equatorial2vec(RA, DE)
        in_circle = inter_star_angle_vec(center_vec, s_vec) < p.theta_FOV
        ID_in_circle_list = list(I[in_circle])
        # select stars randomly
        obs_setID = rand.random_select(
            p.n_obs, ID_in_circle_list, seed=seeds_one[1])
        s_array = equatorial2vec(RA[obs_setID], DE[obs_setID])
        s_list = list(s_array)
        # add noise on star position
        s_hat_list, theta_error_list = add_observation_noise(
            s_list, p, seed=seeds_one[2])
        ### searching graph ###
        candi_setIDs = matching_subgraph(
            p.n_obs, s_hat_list, p.k*p.theta_img, pairstar_db)
        ### calculation stats ###
        # weight
        theta_errors = np.array(theta_error_list)
        weights = calc_weight(theta_errors, p)
        weight = np.prod(weights)
        # stats
        matching_num, multiple, unique, noexist, included = calc_stats(
            candi_setIDs, obs_setID)
        ### LOGGING ###
        logger_case_value.append(
            [matching_num, multiple, unique, noexist, included, weight])
        ### RENDERING ###
        logs = logger_case_value.get_array()
        mean = logs[:, -1].mean()
        print(f'{sample_i+1}/{p.sample_N} : E[weight] = {mean}')
    # save
    logger_case_value.save(
        f'stats_{p.n_obs}_{p.seed}_{p.Vmax}_{p.theta_FOV}_{p.theta_img}_{p.k}_{sampling_type}')


def calc_stats(candi_setIDs, obs_setID):
    # matching number
    matching_num = len(candi_setIDs)
    multiple = matching_num > 1
    unique = matching_num == 1
    noexist = matching_num == 0
    # correct pair is included in candidate pairs
    included = False
    for candi_setID in candi_setIDs:
        if set(obs_setID) == set(candi_setID):
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
        default=2,
        help="ganbare")
    parser.add_argument(
        "--k",
        type=int,
        default=2,
        help="ganbare")
    args = parser.parse_args()
    #
    main(args)
    print('Task complete')
