import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils

# Hyperparameter
sampling_type = 0
sample_N = 10000
seed = 10
n_obs_list = [2, 3, 4, 5, 6]
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
k_list = [2.0**i for i in [-1.0, -0.5, 0.0, 0.5, 1.0]]
# system
log_dir = './log/subgraph_monte/'


def main():
    utils.font_setting()
    #
    ps = list(itertools.product(n_obs_list, Vmax_list,
                                theta_FOV_list, theta_img_list, k_list))
    #
    data_set = {}
    for i, [n_obs, Vmax, theta_FOV, theta_img, k] in enumerate(ps):
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
            pass
        #
        data['time_mean'] = df_obs['time'].mean()
        data['time_75'] = df_obs['time'].quantile(0.75)
        data['time_50'] = df_obs['time'].quantile()
        data['time_25'] = df_obs['time'].quantile(0.25)
        data['time_min'] = df_obs['time'].min()
        data['time_max'] = df_obs['time'].max()
        #
        data['matching_num_mean'] = df_obs['matching_num'].mean()
        data['matching_num_75'] = df_obs['matching_num'].quantile(0.75)
        data['matching_num_50'] = df_obs['matching_num'].quantile()
        data['matching_num_25'] = df_obs['matching_num'].quantile(0.25)
        data['matching_num_min'] = df_obs['matching_num'].min()
        data['matching_num_max'] = df_obs['matching_num'].max()
        #
        data['multiple_num'] = df_obs['multiple'].sum()
        data['unique_num'] = df_obs['unique'].sum()
        data['noexist_num'] = df_obs['noexist'].sum()
        #
        data['included_num'] = df_obs['included'].sum()
        #
        data['collect'] = (df_obs['unique']*df_obs['included']).sum()
        data['miss'] = (df_obs['unique']*(1 - df_obs['included'])).sum()
        data['ambiguous'] = (1 - df_obs['unique']).sum()
        #
        data_set[i] = data
    #
    df = pd.DataFrame.from_dict(data_set, orient="index")
    df.to_csv(log_dir + 'summary.csv')


if __name__ == '__main__':
    main()
    print('task completed')
