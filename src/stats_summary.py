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
    utils.font_setting()
    #
    ps = list(itertools.product(n_obs_list, Vmax_list,
                                theta_FOV_list, theta_img_list, k_list))
    #
    data_set = {}
    for i, [n_obs, Vmax, theta_FOV, theta_img, k] in enumerate(tqdm.tqdm(ps)):
        data = {}
        #
        data['n_obs'] = n_obs
        data['Vmax'] = Vmax
        data['theta_FOV'] = theta_FOV
        data['theta_img'] = theta_img
        data['k'] = k
        data['epsilon'] = theta_img*k
        #
        fname = f'stats_{sample_N}_{seed}_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}_{sampling_type}'
        df = pd.read_csv(log_dir + fname + '.csv', index_col=0)

        ### calc stats ###
        #
        data['calc_num'] = len(df)
        #
        obs_flag = df['observable'] == 1
        data['obs_num'] = obs_flag.sum()
        df_obs = df[obs_flag]
        if len(df_obs) < 1:
            data['time_mean'] = df_obs['time'].mean()
            data['obs_prob'] = data['obs_num'] / data['calc_num']
            #
            data['determined_num'] = None
            data['determined_prob'] = None
            data['correct_num'] = None
            data['correct_prob'] = None
            data['correct_determined_prob'] = None
            data['correct_prob_on_determined'] = None
            data['incorrect_determined_prob'] = None
            data['incorrect_prob_on_determined'] = None
        else:
            data['time_mean'] = df_obs['time'].mean()
            data['obs_prob'] = data['obs_num'] / data['calc_num']
            #
            data['determined_num'] = df_obs['determined'].sum()
            data['determined_prob'] = data['determined_num'] / data['obs_num']
            data['correct_num'] = df_obs['correct'].sum()
            data['correct_prob'] = data['correct_num'] / data['obs_num']
            #
            data['correct_determined_num'] = (
                df_obs['correct']*df_obs['determined']).sum()
            data['correct_determined_prob'] = data['correct_determined_num'] / \
                data['obs_num']
            data['incorrect_determined_num'] = (
                (1 - df_obs['correct'])*df_obs['determined']).sum()
            data['incorrect_determined_prob'] = data['incorrect_determined_num'] / \
                data['obs_num']
            #
            if data['determined_prob'] == 0.0:
                data['correct_prob_on_determined'] = None
                data['incorrect_prob_on_determined'] = None
            else:
                data['correct_prob_on_determined'] = data['correct_determined_prob'] / \
                    data['determined_prob']
                data['incorrect_prob_on_determined'] = data['incorrect_determined_prob'] / \
                    data['determined_prob']
        #
        data_set[i] = data
    #
    df = pd.DataFrame.from_dict(data_set, orient="index")
    df.to_csv(log_dir + f'summary_{sample_N}_{seed}.csv')


if __name__ == '__main__':
    main()
    print('task completed')
