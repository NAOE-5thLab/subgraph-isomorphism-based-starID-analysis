import os
import itertools
import random
import copy

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

import utils

utils.font_setting()


# Hyperparameter
sampling_type = 0
sample_N = 50000
seed = 10
n_obs_list = [2, 3, 4, 5, 6]
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.0, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0]]
k_list = [2.0**i for i in [0.0, 0.5, 1.0]]
# param_col, 'Vmax', 'theta_FOV', 'theta_img', 'k'

# system
log_dir = './log/subgraph_monte/'


def split_df(df, condi_col, condi_value_list):
    df_list = []
    condition_params_list = list(itertools.product(*condi_value_list))
    for condition_params in condition_params_list:
        c = True
        for j in range(len(condi_col)):
            c_temp_low = (df[condi_col[j]] > condition_params[j] -
                          np.abs(condition_params[j])*1.0e-1)
            c_temp_high = (df[condi_col[j]] < condition_params[j] +
                           np.abs(condition_params[j])*1.0e-1)
            c = c & c_temp_low & c_temp_high
        df_list.append(df[c])
    return condition_params_list, df_list


def plot_tar(df_summary, target_col, condi_col, condi_value_list, tag=''):
    scale = 3
    fig = plt.figure(figsize=(scale*len(condi_col), scale))
    axes = []
    for i in range(len(condi_col)):
        axes.append(fig.add_subplot(1, len(condi_col), i+1))

    for i, param_col in enumerate(condi_col):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_summary, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(df[param_col], df[target_col], ls='-', lw=0.2)
            if param_col == 'theta_img':
                axes[i].set_xscale('log')
            axes[i].set_ylabel(target_col)
            axes[i].set_xlabel(param_col)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/{tag}'):
        os.makedirs(f'./img/{tag}')
    fig.savefig(f'./img/{tag}/{target_col}.pdf')
    plt.close(fig)


def plot_hist(df_summary, condi_col, tag=''):
    scale = 3
    fig = plt.figure(figsize=(scale*len(condi_col), scale))
    axes = []
    for i in range(len(condi_col)):
        axes.append(fig.add_subplot(1, len(condi_col), i+1))

    for i, param_col in enumerate(condi_col):
        if param_col == 'theta_img':
            axes[i].hist(df_summary[param_col], bins=np.logspace(
                np.log10(1.0e-7), np.log10(2.0e-3), 20))
            axes[i].set_xscale('log')
        else:
            axes[i].hist(df_summary[param_col])
        if param_col == 'Vmax':
            axes[i].set_xlim(1.0, 6.0)
        axes[i].set_xlabel(param_col)

    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/{tag}'):
        os.makedirs(f'./img/{tag}')
    fig.savefig(f'./img/{tag}/hist.pdf')
    plt.close(fig)


def main():
    # Load
    df_summary = pd.read_csv(
        log_dir + f'summary_{sample_N}_{seed}.csv', index_col=0)
    df_summary.head()    
    # 
    for (Vmax, theta_FOV) in list(itertools.product(Vmax_list, theta_FOV_list)):
        # condition
        tag = f'Vmax{Vmax}theta_FOV{theta_FOV}'
        # 
        vmax_flag_low = (df_summary['Vmax'] > Vmax - np.abs(Vmax)*1.0e-1)
        vmax_flag_high = (df_summary['Vmax'] < Vmax + np.abs(Vmax)*1.0e-1)
        theta_FOV_flag_low = (df_summary['theta_FOV'] > theta_FOV - np.abs(theta_FOV)*1.0e-1)
        theta_FOV_flag_high = (df_summary['theta_FOV'] < theta_FOV + np.abs(theta_FOV)*1.0e-1)
        df = df_summary[
            (vmax_flag_low&vmax_flag_high)&(theta_FOV_flag_low&theta_FOV_flag_high)]    
        # plot
        condi_col = ['n_obs', 'theta_img', 'k']
        condi_value_list = [n_obs_list, theta_img_list, k_list]
        #
        plot_hist(df, condi_col, tag=tag)
        #
        plot_tar(df, 'time_mean', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'matching_num_mean', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'obs_prob', condi_col, condi_value_list, tag=tag)
        # 
        plot_tar(df, 'unique_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'ambiguous_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'multiple_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'noexist_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'included_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'notincluded_prob', condi_col, condi_value_list, tag=tag)
        # 
        plot_tar(df, 'included_unique_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'notincluded_unique_num', condi_col, condi_value_list, tag=tag)
        # 
        plot_tar(df, 'included_on_unique_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'notincluded_on_unique_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'unique_on_included_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'multiple_on_included_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'noexist_on_included_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'unique_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'multiple_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
        plot_tar(df, 'noexist_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
    

if __name__ == '__main__':
    main()
    print('task completed')
