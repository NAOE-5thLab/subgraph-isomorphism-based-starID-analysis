import os
import itertools
import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils

# Hyperparameter
sampling_type = 0
sample_N = 10000
seed = 10
n_obs_list = [2, 3, 4, 5, 6, 7, 8]
# DB
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
# simulation
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
# subgraph matching
k_list = [2.0**i for i in [-0.5, 0.0, 0.5, 1.0]]
# system
log_dir = './log/subgraph_monte/'


def main():
    #
    ps = list(itertools.product(
        range(len(Vmax_list)),
        range(len(theta_FOV_list)),
        range(len(theta_img_list)),
        range(len(k_list))
    ))
    #
    exist = []
    not_exist = []
    for i, indices in enumerate(tqdm.tqdm(ps)):
        Vmax = Vmax_list[indices[0]]
        theta_FOV = theta_FOV_list[indices[1]]
        theta_img = theta_img_list[indices[2]]
        k = k_list[indices[3]]
        #
        fname = f'stats_{sample_N}_{seed}_{n_obs_list[-1]}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}_{sampling_type}'
        if os.path.isfile(log_dir + fname + '.csv'):
            exist.append(indices)
        else:
            not_exist.append(indices)
    #
    print(exist)
    print(not_exist)


if __name__ == '__main__':
    main()
    print('task completed')
