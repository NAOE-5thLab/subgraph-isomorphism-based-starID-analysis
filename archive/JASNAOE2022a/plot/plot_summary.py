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


def filter_df(df_summary):
    condi_cols = ['n_obs', 'Vmax', 'theta_img', 'k']
    condi_value_lists = [n_obs_list, Vmax_list, theta_img_list, k_list]
    condition_params_list, df_list = split_df(
        df_summary, condi_cols, condi_value_lists)
    #
    params_list = []
    for i, df in enumerate(df_list):
        corr = df[['theta_FOV', 'multiple_prob']].corr().iloc[0, 1]
        if corr < -0.5:
            params_list.append(condition_params_list[i])
    #
    f = False
    for params in params_list:
        c = True
        for i in range(len(condi_cols)):
            c_temp_low = (df_summary[condi_cols[i]] > params[i] -
                          np.abs(params[i])*1.0e-1)
            c_temp_high = (df_summary[condi_cols[i]] < params[i] +
                           np.abs(params[i])*1.0e-1)
            c = c & c_temp_low & c_temp_high
        f = f | c
    return df_summary[f]  # [f] [~f]


def main():
    # Load
    df_summary = pd.read_csv(
        log_dir + f'summary_{sample_N}_{seed}.csv', index_col=0)
    df_summary.head()
    # condition
    # tag = ''
    # tag = 'Vmag1.5'
    # df_summary = df_summary[df_summary['Vmax'] == 1.5]
    # tag = 'Vmag3.5'
    # df_summary = df_summary[df_summary['Vmax'] == 3.5]
    # tag = 'k1'
    # df_summary = df_summary[df_summary['k'] == 1]
    # tag = 'Vmag1.5k1'
    # df_summary = df_summary[(df_summary['Vmax'] == 1.5) & (df_summary['k'] == 1)]
    tag = 'FOVcorr<-0.5'
    df_summary = filter_df(df_summary)

    # plot
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'theta_img', 'k']
    condi_value_list = [n_obs_list, Vmax_list, theta_FOV_list, theta_img_list, k_list]
    #
    plot_hist(df_summary, condi_col, tag=tag)
    #
    plot_tar(df_summary, 'time_mean', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'matching_num_mean', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'obs_prob', condi_col, condi_value_list, tag=tag)
    # 
    plot_tar(df_summary, 'unique_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'ambiguous_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'multiple_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'noexist_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'included_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'notincluded_prob', condi_col, condi_value_list, tag=tag)
    # 
    plot_tar(df_summary, 'included_unique_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'notincluded_unique_num', condi_col, condi_value_list, tag=tag)
    # 
    plot_tar(df_summary, 'included_on_unique_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'notincluded_on_unique_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'unique_on_included_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'multiple_on_included_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'noexist_on_included_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'unique_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'multiple_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
    plot_tar(df_summary, 'noexist_on_notincluded_prob', condi_col, condi_value_list, tag=tag)
    

if __name__ == '__main__':
    main()
    print('task completed')
