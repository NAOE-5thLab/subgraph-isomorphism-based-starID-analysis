import os
import itertools
import random
import copy

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from config import Param
import utils

utils.font_setting()


conf = Param()
epsilon_list = [pair[0]*pair[1]
            for pair in list(itertools.product(conf.theta_img_list, conf.k_list))]
epsilon_list = sorted(list(set(epsilon_list)))


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


def plot_time(df_summary):
    scale = 1.5
    fig = plt.figure(figsize=(scale*1.6*2, scale*2))
    axes = [fig.add_subplot(2, 2, i+1) for i in range(4)]

    target_col = 'time_mean'
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'epsilon']
    condi_value_list = [conf.n_obs_list, conf.Vmax_list, conf.theta_FOV_list, epsilon_list]
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
    condi_value_list = [conf.n_obs_list, conf.Vmax_list, conf.theta_FOV_list, epsilon_list]
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
        conf.n_obs_list, conf.Vmax_list, conf.theta_FOV_list, conf.theta_img_list, conf.k_list]
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

    tag = 'k1andn_obs2'
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

    tag = 'k1orn_obs2'
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
        conf.n_obs_list, conf.Vmax_list, conf.theta_FOV_list, conf.theta_img_list, conf.k_list]
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
        conf.n_obs_list, conf.Vmax_list, conf.theta_FOV_list, conf.theta_img_list, conf.k_list]
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
        conf.log_dir + f'summary.csv', index_col=0)
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
