import numpy as np

import rand


class Param:
    ### Hyperparameter ###
    # DB
    Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
    theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
    # simulation
    theta_img_list = [np.pi/180*10**i for i in [
        -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
    # subgraph matching
    k_list = [2.0**i for i in [-0.5, 0.0, 0.5, 1.0]]
    n_obs_list = [2, 3, 4, 5, 6, 7, 8]

    ### basic parameter ###
    # system
    log_dir = './log/subgraph_monte/'
    log_temp_dir = './log/subgraph_monte/temp/'
    # simulation
    calc_seeds = range(10000)
    parallel_num = 10
    alpha = 0.99
