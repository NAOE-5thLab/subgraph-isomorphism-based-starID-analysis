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


n_obs_list = [2, 3, 4, 5, 6]
# DB
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
# simulation
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.0, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0]]
# subgraph matching
k_list = [2.0**i for i in [0.0, 0.5, 1.0]]

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


def plot_obs_prob(df_summary, tag=''):
    scale = 3
    fig = plt.figure(figsize=(scale*1.6, scale))

    df = df_summary
    axes = [
        fig.add_subplot(2, 2, 1), fig.add_subplot(2, 2, 2),
        fig.add_subplot(2, 2, 3), fig.add_subplot(2, 2, 4)
    ]
    ls_list = ["-", "--", ":", "-.", (0, (3, 1, 1, 1, 1, 1))]
    marker_list = ['o', '^', 's', 'd', 'x']
    # color_list = ['#267E63', '#38AEBC', '#4B9EFB', '#787FFC', '#C5A5FD']

    # plot
    for j, theta_FOV in enumerate(theta_FOV_list):
        c_low = (df['theta_FOV'] > theta_FOV - np.abs(theta_FOV)*1.0e-1)
        c_high = (df['theta_FOV'] < theta_FOV + np.abs(theta_FOV)*1.0e-1)
        df_1 = df[c_low & c_high]
        #
        for i, n_obs in enumerate(n_obs_list[::2]):
            c_low = (df_1['n_obs'] > n_obs - np.abs(n_obs)*1.0e-1)
            c_high = (df_1['n_obs'] < n_obs + np.abs(n_obs)*1.0e-1)
            df_2 = df_1[c_low & c_high]
            axes[j].plot(
                df_2['Vmax'], df_2['obs_prob'],
                color='black', ls=ls_list[i], lw=0.7,
                marker=marker_list[i], markersize=2,
                label='$p = $'+f'{n_obs}')
        axes[j].set_xlim(1.4, 5.6)
        axes[j].set_ylim(-0.05, 1.05)
        axes[j].set_title(
            '$\\theta_{\mathrm{FOV}} = $'+f'{int(theta_FOV*180/np.pi+0.5)}', fontsize=9)
        # axes[j].text(1.5, 1.1, '$\\theta_{\mathrm{FOV}} = $'+f'{int(theta_FOV*180/np.pi+0.5)}', fontsize=9)
        if j > 1:
            axes[j].set_xlabel('$V_{\mathrm{max}}$')
        else:
            axes[j].set_xticklabels([])
        if j % 2 == 0:
            axes[j].set_ylabel('obs_prob')
        else:
            axes[j].set_yticklabels([])
        if j == 0:
            axes[j].legend(fontsize='x-small')
    #
    fig.tight_layout()
    fig.savefig(f'./img/{tag}/prop_obs.pdf')
    plt.close(fig)


def plot_time_and_matching_num(df_summary, tag=''):
    scale = 3
    fig = plt.figure(figsize=(scale*1.6*2, scale))
    axes = [
        [
            fig.add_subplot(2, 4, 1), fig.add_subplot(2, 4, 2),
            fig.add_subplot(2, 4, 3), fig.add_subplot(2, 4, 4),
        ], [
            fig.add_subplot(2, 4, 5), fig.add_subplot(2, 4, 6),
            fig.add_subplot(2, 4, 7), fig.add_subplot(2, 4, 8)
        ]
    ]
    #
    condi_col = ['n_obs', 'Vmax', 'theta_FOV', 'epsilon']
    condi_value_list = [n_obs_list, Vmax_list, theta_FOV_list, epsilon_list]
    tar_col = ['time_mean', 'matching_num_mean']

    xlabels = ['$p$', '$V_{\mathrm{max}}$',
               '$\\theta_{\mathrm{FOV}}$ [deg.]', '$\\varepsilon$ [deg.]']
    ylabels = ['$\\overline{T}$', '$\\overline{|C_{12\\cdots p}|}$']

    for j in range(2):
        for i in range(4):
            max_data = []
            min_data = []
            for value in condi_value_list[i]:
                c_low = (df_summary[condi_col[i]] >
                         value - np.abs(value)*1.0e-1)
                c_high = (df_summary[condi_col[i]] <
                          value + np.abs(value)*1.0e-1)
                df_1 = df_summary[c_low & c_high]
                #
                max_data.append(df_1[tar_col[j]].max())
                min_data.append(df_1[tar_col[j]].min())
            max_data = np.array(max_data)
            min_data = np.array(min_data)
            if i >= 2:
                x = np.array(condi_value_list[i])*180/np.pi
            else:
                x = np.array(condi_value_list[i])
            #
            axes[j][i].plot(x, max_data, ls='-', lw=0.7, marker='o',
                            markersize=2, color='black', label='max')
            axes[j][i].plot(x, min_data, ls='--', lw=0.7, marker='^',
                            markersize=2, color='black', label='min')
            if i == 3:
                axes[j][i].set_xscale('log')
            if i % 4 == 0:
                axes[j][i].set_ylabel(ylabels[j])
            else:
                axes[j][i].set_yticklabels([])
            if j == 1:
                axes[j][i].set_xlabel(xlabels[i])
            if i == 0 & j == 0:
                axes[j][i].legend(fontsize='x-small')
    #

    #
    fig.tight_layout()
    fig.savefig(f'./img/{tag}/time_and_matching_num.pdf')


def main():
    # Load
    df_summary = pd.read_csv(
        log_dir + f'summary_{sample_N}_{seed}.csv', index_col=0)
    df_summary.head()
    # condition
    tag = ''
    #
    plot_obs_prob(df_summary, tag=tag)
    plot_time_and_matching_num(df_summary, tag=tag)


if __name__ == '__main__':
    main()
    print('task completed')

    # plot and legend
    # ls_list = ["-", "--", ":", "-."]
    # marker_list = ['o', '^', 'x', 'd', 's']
    # color_list = ['#267E63', '#38AEBC', '#4B9EFB', '#787FFC', '#C5A5FD']
    # for i, n_obs in enumerate(n_obs_list):
    #     c_low = (df['n_obs'] > n_obs - np.abs(n_obs)*1.0e-1)
    #     c_high = (df['n_obs'] < n_obs + np.abs(n_obs)*1.0e-1)
    #     df_1 = df[c_low&c_high]
    #     for j, theta_FOV in enumerate(theta_FOV_list):
    #         c_low = (df_1['theta_FOV'] > theta_FOV - np.abs(theta_FOV)*1.0e-1)
    #         c_high = (df_1['theta_FOV'] < theta_FOV + np.abs(theta_FOV)*1.0e-1)
    #         df_2 = df_1[c_low&c_high]
    #         ax.plot(
    #             df_2['Vmax'], df_2['obs_prob'], color = color_list[j],
    #             ls=ls_list[j], lw=1)
    #         ax.plot(
    #             df_2['Vmax'], df_2['obs_prob'], color = 'Black',
    #             ls='None', marker=marker_list[i], markersize=2)
    # ax.set_xlim(1,6)
    # ax.set_ylim(-0.05, 1.05)
    # ax.set_xlabel('$V_{\mathrm{max}}$')
    # ax.set_ylabel('obs_prob')
    # for i, n_obs in enumerate(n_obs_list):
    #     ax.plot(
    #         -1, -1, color = 'Black',
    #         ls='None', marker=marker_list[i], markersize=2,
    #         label='$n_{\mathrm{obs}} = $'+f'{n_obs}')
    # for j, theta_FOV in enumerate(theta_FOV_list):
    #     ax.plot(
    #         df_2['Vmax'], df_2['obs_prob'], color = color_list[j],
    #         ls=ls_list[j], lw=1,
    #         label='$\\theta_{\mathrm{FOV}} = $'+f'{int(theta_FOV*180/np.pi+0.5)}')
    # ax.legend()

    # scatter and colorbar
    # cm = plt.cm.get_cmap('RdYlBu')
    # sca = ax.scatter(df['Vmax'], df['theta_FOV'], c=df['obs_prob'], vmin=0, vmax=1, cmap=cm)
    # ax.set_xlabel('$V_{\mathrm{max}}$')
    # ax.set_ylabel('$\\theta_{\mathrm{FOV}}$')
    # fig.colorbar(sca, ax=ax)
