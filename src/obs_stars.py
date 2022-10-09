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

    ### basic parameter ###
    # system
    log_dir = './log/'
    seed = 10
    # simulation
    sample_N = 10000

    def __init__(self, args):
        ### Hyperparameter ###
        self.Vmax = self.Vmax_list[args.Vmax]
        self.theta_FOV = self.theta_FOV_list[args.theta_FOV]


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

    ### Prepair LOGGER ###
    header = ['incircle']
    logger_case_value = log.Logger(p.log_dir, header)
    logger_case_value.reset()

    ### MONTE CARLO SIMULATION ###
    # seed
    seeds_monte = gen_seeds(p.seed, p.sample_N)
    for sample_i in range(p.sample_N):
        # select center randomly
        center_vec = rand.uniform_spherical_vector(seed=seeds_monte[sample_i])
        # collect stars in circle (FOV)
        s_vec = equatorial2vec(RA, DE)
        in_circle = inter_star_angle_vec(center_vec, s_vec) < p.theta_FOV
        ID_in_circle_list = list(I[in_circle])
        # num
        in_circle_num = len(ID_in_circle_list)

        ### LOGGING ###
        logger_case_value.append([in_circle_num])
    # save
    logger_case_value.save(f'obs_stars_{p.Vmax}_{p.theta_FOV}')


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
    args = parser.parse_args()
    #
    main(args)
    print('Task complete')
