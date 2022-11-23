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
sample_N = 10000
seed = 10
# n_obs_list = [2, 3, 4, 5, 6, 7, 8]
# # DB
# Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
# theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
# # simulation
# theta_img_list = [np.pi/180*10**i for i in [
#     -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
# # subgraph matching
# k_list = [2.0**i for i in [-0.5, 0.0, 0.5, 1.0]]
# param_col, 'Vmax', 'theta_FOV', 'theta_img', 'k'


n_obs_list = [2, 3, 4, 5, 6, 7, 8]
# DB
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
# simulation
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
# subgraph matching
k_list = [2.0**i for i in [-0.5, 0.0, 0.5, 1.0]]

epsilon_list = [pair[0]*pair[1]
                for pair in list(itertools.product(theta_img_list, k_list))]
epsilon_list = sorted(list(set(epsilon_list)))

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
            axes[i].plot(df[param_col], df[target_col],
                         color='black', ls='-', lw=0.2)
            # axes[i].plot(df[param_col], df[target_col], ls='-', lw=0.2)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            axes[i].set_ylabel(target_col)
            axes[i].set_xlabel(param_col)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/{tag}'):
        os.makedirs(f'./img/{tag}')
    fig.savefig(f'./img/{tag}/{target_col}.pdf')
    plt.close(fig)


def plot_tar_2(
        df_summary, target_col,
        condi_col, condi_value_list,
        split_col, split_value, split_color, tag=''):
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
        for j, value in enumerate(split_value):
            c_temp_low = (df_summary[split_col] > value - np.abs(value)*1.0e-1)
            c_temp_high = (df_summary[split_col] <
                           value + np.abs(value)*1.0e-1)
            df_split = df_summary[c_temp_low & c_temp_high]
            #
            condition_params_list, df_list = split_df(
                df_split, condi_cols, condi_value_lists)
            for df in df_list:
                axes[i].plot(df[param_col], df[target_col],
                             color=split_color[j], ls='-', lw=0.2)
                if param_col == 'theta_img' or param_col == 'epsilon':
                    axes[i].set_xscale('log')
                axes[i].set_ylabel(target_col)
                axes[i].set_xlabel(param_col)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/{tag}'):
        os.makedirs(f'./img/{tag}')
    fig.savefig(f'./img/{tag}/{target_col}.pdf')
    plt.close(fig)


def plot_tar_3(df_summary, target_col, target_col_2, condi_col, condi_value_list, tag=''):
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
            axes[i].plot(df[target_col_2], df[target_col],
                         color='black', ls='-', lw=0.2)
            # axes[i].plot(df[param_col], df[target_col], ls='-', lw=0.2)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            axes[i].set_title(param_col)
            axes[i].set_ylabel(target_col)
            axes[i].set_xlabel(target_col_2)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/{tag}'):
        os.makedirs(f'./img/{tag}')
    fig.savefig(f'./img/{tag}/{target_col}_{target_col_2}.pdf')
    plt.close(fig)


def plot_time(df_summary):
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*2))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(4)]

    target_col = 'time_mean'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'epsilon']
    condi_value_list = [n_obs_list, Vmax_list, theta_FOV_list, epsilon_list]
    labels = ['$p$', '$V_{\mathrm{max}}$',
              '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\varepsilon$']

    for i, param_col in enumerate(condi_col):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_summary, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$T_{\mathrm{cal}}$ [s]')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.5, 13.0)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}.pdf')
    plt.close(fig)


def plot_obs_prob(df_summary):
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*2))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(4)]

    target_col = 'obs_prob'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'epsilon']
    condi_value_list = [n_obs_list, Vmax_list, theta_FOV_list, epsilon_list]
    labels = ['$p$', '$V_{\mathrm{max}}$',
              '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\varepsilon$']

    for i, param_col in enumerate(condi_col):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_summary, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$P(O)$')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.05, 1.05)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}.pdf')
    plt.close(fig)


def plot_determined_prob(df_summary):
    target_col = 'determined_prob'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'theta_img', 'k']
    condi_value_list = [
        n_obs_list, Vmax_list, theta_FOV_list, theta_img_list, k_list]
    labels = ['$p$', '$V_{\mathrm{max}}$',
              '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\theta_{\mathrm{img}}$ [deg.]', '$k$']
    #
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*3))
    axes = [fig.add_subplot(3, 2, i+1) for i in range(5)]

    df_condi = df_summary
    for i, param_col in enumerate(condi_col):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_condi, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$P(A)$')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.05, 1.05)
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}.pdf')
    plt.close(fig)

    ### k < 1 and n_obs > 2 ###
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*2))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(4)]

    tag = 'k<1andn_obs>2'
    flag = (df_summary['k'] < 1.01) & (df_summary['n_obs'] > 2.5)
    df_condi = df_summary[flag]

    for i, param_col in enumerate(condi_col[:4]):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_condi, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$P(A)$')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.05, 1.05)
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}_{tag}.pdf')
    plt.close(fig)

    ### k > 1 or n_obs = 2 ###
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*2))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(4)]

    tag = 'k>1orn_obs=2'
    flag = (df_summary['k'] < 1.01) & (df_summary['n_obs'] > 2.5)
    df_condi = df_summary[~flag]

    for i, param_col in enumerate(condi_col[:4]):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_condi, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$P(A)$')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.05, 1.05)
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}_{tag}.pdf')
    plt.close(fig)


def plot_correct_prob_on_determined(df_summary):
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*3))
    axes = [fig.add_subplot(3, 2, i+1) for i in range(5)]

    target_col = 'correct_prob_on_determined'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'theta_img', 'k']
    condi_value_list = [
        n_obs_list, Vmax_list, theta_FOV_list, theta_img_list, k_list]
    labels = ['$p$', '$V_{\mathrm{max}}$',
              '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\theta_{\mathrm{img}}$ [deg.]', '$k$']

    for i, param_col in enumerate(condi_col):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_summary, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[param_col], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            if i % 2 == 0:
                axes[i].set_ylabel('$P_{A}(B)$')
            axes[i].set_xlabel(labels[i])
            axes[i].set_ylim(-0.05, 1.05)
    #
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}.pdf')
    plt.close(fig)


def pair_plot_determined_prob_obs_prob(df_summary):
    target_col = 'determined_prob'
    target_col_2 = 'obs_prob'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'theta_img', 'k']
    condi_value_list = [
        n_obs_list, Vmax_list, theta_FOV_list, theta_img_list, k_list]
    labels = ['$p$', '$V_{\mathrm{max}}$',
              '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\theta_{\mathrm{img}}$ [deg.]', '$k$']
    #
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*3))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(3)]

    flag = (df_summary['k'] < 1.01) & (df_summary['n_obs'] > 2.5)
    df_condi = df_summary[~flag]
    for i, param_col in enumerate(condi_col[:3]):
        condi_cols = copy.copy(condi_col)
        condi_cols.pop(i)
        condi_value_lists = copy.copy(condi_value_list)
        condi_value_lists.pop(i)
        #
        condition_params_list, df_list = split_df(
            df_condi, condi_cols, condi_value_lists)
        for df in df_list:
            axes[i].plot(
                df[target_col_2], df[target_col],
                color='black', ls='-', alpha=0.5, lw=0.3,
                marker='.', markersize=1)
            if param_col == 'theta_img' or param_col == 'epsilon':
                axes[i].set_xscale('log')
            axes[i].set_ylabel('$P(A)$')
            axes[i].set_xlabel('$P(O)$')
            axes[i].set_title(labels[i])
            axes[i].set_aspect('equal')
    fig.tight_layout()
    if not os.path.isdir(f'./img/paper'):
        os.makedirs(f'./img/paper')
    fig.savefig(f'./img/paper/{target_col}_{target_col_2}.pdf')
    plt.close(fig)


def main():
    # Load
    df_summary = pd.read_csv(
        log_dir + f'summary_{sample_N}_{seed}.csv', index_col=0)
    df_summary.head()
    #
    plot_time(df_summary)
    plot_obs_prob(df_summary)
    plot_determined_prob(df_summary)
    plot_correct_prob_on_determined(df_summary)
    pair_plot_determined_prob_obs_prob(df_summary)


if __name__ == '__main__':
    main()
    print('task completed')
