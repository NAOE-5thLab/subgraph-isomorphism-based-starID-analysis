import itertools
import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils

# Hyperparameter
sampling_type = 0
sample_N = 50000
seed = 10
n_obs_list = [2, 3, 4, 5, 6]
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
            data['matching_num_mean'] = df_obs['matching_num'].mean()
            #
            data['obs_prob'] = data['obs_num'] / data['calc_num']
            #
            data['unique_num'] = None
            data['unique_prob'] = None
            data['ambiguous_num'] = None
            data['ambiguous_prob'] = None
            data['multiple_num'] = None
            data['multiple_prob'] = None
            data['noexist_num'] = None
            data['noexist_prob'] = None
            #
            data['included_num'] = None
            data['included_prob'] = None
            data['notincluded_num'] = None
            data['notincluded_prob'] = None
            #
            data['included_unique_num'] = None
            data['included_unique_prob'] = None
            data['notincluded_unique_num'] = None
            data['notincluded_unique_prob'] = None
            data['notincluded_noexist_num'] = None
            data['notincluded_noexist_prob'] = None
            #
            data['included_on_unique_prob'] = None
            data['notincluded_on_unique_prob'] = None
            data['noexist_on_notincluded_prob'] = None
        else:
            data['time_mean'] = df_obs['time'].mean()
            data['matching_num_mean'] = df_obs['matching_num'].mean()
            #
            data['obs_prob'] = data['obs_num'] / data['calc_num']
            #
            data['unique_num'] = df_obs['unique'].sum()
            data['unique_prob'] = data['unique_num'] / data['obs_num']
            data['ambiguous_num'] = (1 - df_obs['unique']).sum()
            data['ambiguous_prob'] = data['ambiguous_num'] / data['obs_num']
            data['multiple_num'] = df_obs['multiple'].sum()
            data['multiple_prob'] = data['multiple_num'] / data['obs_num']
            data['noexist_num'] = df_obs['noexist'].sum()
            data['noexist_prob'] = data['noexist_num'] / data['obs_num']
            #
            data['included_num'] = df_obs['included'].sum()
            data['included_prob'] = data['included_num'] / data['obs_num']
            data['notincluded_num'] = (1-df_obs['included']).sum()
            data['notincluded_prob'] = data['notincluded_num'] / data['obs_num']
            #
            data['included_unique_num'] = (df_obs['unique']*df_obs['included']).sum()
            data['included_unique_prob'] = data['included_unique_num'] / data['obs_num']
            data['included_multiple_num'] = (df_obs['multiple']*df_obs['included']).sum()
            data['included_multiple_prob'] = data['included_multiple_num'] / data['obs_num']
            data['included_noexist_num'] = (df_obs['noexist']*df_obs['included']).sum()
            data['included_noexist_prob'] = data['included_noexist_num'] / data['obs_num']
            data['notincluded_unique_num'] = (df_obs['unique']*(1 - df_obs['included'])).sum()
            data['notincluded_unique_prob'] = data['notincluded_unique_num'] / data['obs_num']
            data['notincluded_multiple_num'] = (df_obs['multiple']*(1 - df_obs['included'])).sum()
            data['notincluded_multiple_prob'] = data['notincluded_multiple_num'] / data['obs_num']
            data['notincluded_noexist_num'] = (df_obs['noexist']*(1 - df_obs['included'])).sum()
            data['notincluded_noexist_prob'] = data['notincluded_noexist_num'] / data['obs_num']
            # 
            data['included_on_unique_prob'] = data['included_unique_prob'] / data['unique_prob']
            data['notincluded_on_unique_prob'] = data['notincluded_unique_prob'] / data['unique_prob']
            #
            data['unique_on_included_prob'] = data['included_unique_prob'] / data['included_prob']
            data['multiple_on_included_prob'] = data['included_multiple_prob'] / data['included_prob']
            data['noexist_on_included_prob'] = data['included_noexist_prob'] / data['included_prob']
            data['unique_on_notincluded_prob'] = data['notincluded_unique_prob'] / data['notincluded_prob']
            data['multiple_on_notincluded_prob'] = data['notincluded_multiple_prob'] / data['notincluded_prob']
            data['noexist_on_notincluded_prob'] = data['notincluded_noexist_prob'] / data['notincluded_prob']
        #
        data_set[i] = data
    #
    df = pd.DataFrame.from_dict(data_set, orient="index")
    df.to_csv(log_dir + f'summary_{sample_N}_{seed}.csv')


if __name__ == '__main__':
    main()
    print('task completed')
