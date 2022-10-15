import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils

utils.font_setting()

# Hyperparameter
sampling_type = 1
sample_N = 100
seed = 10
n_obs_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5]]
k_list = [2.0**i for i in [-1.0, -0.5, 0.0, 0.5, 1.0]]
# system
log_dir = './log/subgraph_monte/'


N_n_obs_list = len(n_obs_list)
N_Vmax_list = len(Vmax_list)
N_theta_FOV_list = len(theta_FOV_list)
N_theta_img_list = len(theta_img_list)
N_k_list = len(k_list)

params = list(
    itertools.product(
        range(N_n_obs_list), range(N_Vmax_list), range(N_theta_FOV_list),
        range(N_theta_img_list), range(N_k_list)))
time_mean = np.zeros(
    [N_n_obs_list, N_Vmax_list, N_theta_FOV_list, N_theta_img_list, N_k_list])
matching_num_mean = np.zeros(
    [N_n_obs_list, N_Vmax_list, N_theta_FOV_list, N_theta_img_list, N_k_list])
unique_mean = np.zeros(
    [N_n_obs_list, N_Vmax_list, N_theta_FOV_list, N_theta_img_list, N_k_list])
included_mean = np.zeros(
    [N_n_obs_list, N_Vmax_list, N_theta_FOV_list, N_theta_img_list, N_k_list])


def process(df):
    observable = df['observable'].to_numpy()
    time = df['time'].to_numpy()
    matching_num = df['matching_num'].to_numpy()
    multiple = df['multiple'].to_numpy()
    unique = df['unique'].to_numpy()
    noexist = df['noexist'].to_numpy()
    included = df['included'].to_numpy()
    weight = df['weight'].to_numpy()
    #
    obs_flag = observable == 1
    overtime_flag = len(observable) < sample_N
    #
    prob_observe = observable.mean()

    multiple[multiple == 'True'] = 1
    multiple[multiple == 'False'] = 0


for param in params:
    #
    n_obs = n_obs_list[param[0]]
    Vmax = Vmax_list[param[1]]
    theta_FOV = theta_FOV_list[param[2]]
    theta_img = theta_img_list[param[3]]
    k = k_list[param[4]]
    #
    fname = f'stats_{sample_N}_{seed}_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}_{sampling_type}'
    df = pd.read_csv(log_dir + fname + '.csv', index_col=0)
    process(df)
    #
    observable = df['observable'].to_numpy()
    time = df['time'].to_numpy()
    matching_num = df['matching_num'].to_numpy()
    multiple = df['multiple'].to_numpy()
    unique = df['unique'].to_numpy()
    noexist = df['noexist'].to_numpy()
    included = df['included'].to_numpy()
    weight = df['weight'].to_numpy()
    #
    obs_flag = observable == 1
    overtime_flag = len(observable) < sample_N
    #
    if overtime_flag:
        time_mean[param[0], param[1], param[2], param[3], param[4]] = 60.0
    else:
        time_mean[param[0], param[1], param[2], param[3], param[4]] = (
            time[obs_flag]*weight[obs_flag]).mean()
    #
    matching_num_mean[param[0], param[1], param[2], param[3], param[4]] = (
        matching_num[obs_flag]*weight[obs_flag]).mean()
    unique_mean[param[0], param[1], param[2], param[3], param[4]] = (
        unique[obs_flag]*weight[obs_flag]).mean()
    included_mean[param[0], param[1], param[2], param[3], param[4]] = (
        included[obs_flag]*weight[obs_flag]).mean()
