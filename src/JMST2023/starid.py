import os
import argparse
import math
import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import special_ortho_group

from database import YaleStarCatalog, StarCatalog, InterAngleCatalog
from matching.subgraph import matching_set_for_analysis
from utils.unit_convert import deg2rad
from utils.seeds import gen_seeds
from utils.spherical_ramdom_vector import limit_uniform_spherical_vector
from utils.rotation import rodrigues_rotation_matrix
from utils.font import font_setting
from utils.inter_star_angle import inter_star_angle_vec

parser = argparse.ArgumentParser()
parser.add_argument(
    "--i",
    type=int,
    default=1,
    help="degree")
args = parser.parse_args()

patterns = [
    {'FOV_deg': 20, 'obsMv': 5.5, 'beta': 0.0},
    {'FOV_deg': 40, 'obsMv': 5.5, 'beta': 0.0},
    {'FOV_deg': 60, 'obsMv': 5.5, 'beta': 0.0},
    {'FOV_deg': 80, 'obsMv': 5.5, 'beta': 0.0},
    {'FOV_deg': 20, 'obsMv': 4.5, 'beta': 0.0},
    {'FOV_deg': 40, 'obsMv': 4.5, 'beta': 0.0},
    {'FOV_deg': 60, 'obsMv': 4.5, 'beta': 0.0},
    {'FOV_deg': 80, 'obsMv': 4.5, 'beta': 0.0},
    {'FOV_deg': 20, 'obsMv': 3.5, 'beta': 0.0},
    {'FOV_deg': 40, 'obsMv': 3.5, 'beta': 0.0},
    {'FOV_deg': 60, 'obsMv': 3.5, 'beta': 0.0},
    {'FOV_deg': 80, 'obsMv': 3.5, 'beta': 0.0},
    {'FOV_deg': 20, 'obsMv': 5.5, 'beta': 0.2},
    {'FOV_deg': 40, 'obsMv': 5.5, 'beta': 0.2},
    {'FOV_deg': 60, 'obsMv': 5.5, 'beta': 0.2},
    {'FOV_deg': 80, 'obsMv': 5.5, 'beta': 0.2},
    {'FOV_deg': 20, 'obsMv': 4.5, 'beta': 0.2},
    {'FOV_deg': 40, 'obsMv': 4.5, 'beta': 0.2},
    {'FOV_deg': 60, 'obsMv': 4.5, 'beta': 0.2},
    {'FOV_deg': 80, 'obsMv': 4.5, 'beta': 0.2},
    {'FOV_deg': 20, 'obsMv': 3.5, 'beta': 0.2},
    {'FOV_deg': 40, 'obsMv': 3.5, 'beta': 0.2},
    {'FOV_deg': 60, 'obsMv': 3.5, 'beta': 0.2},
    {'FOV_deg': 80, 'obsMv': 3.5, 'beta': 0.2},
    {'FOV_deg': 20, 'obsMv': 5.5, 'beta': 0.4},
    {'FOV_deg': 40, 'obsMv': 5.5, 'beta': 0.4},
    {'FOV_deg': 60, 'obsMv': 5.5, 'beta': 0.4},
    {'FOV_deg': 80, 'obsMv': 5.5, 'beta': 0.4},
    {'FOV_deg': 20, 'obsMv': 4.5, 'beta': 0.4},
    {'FOV_deg': 40, 'obsMv': 4.5, 'beta': 0.4},
    {'FOV_deg': 60, 'obsMv': 4.5, 'beta': 0.4},
    {'FOV_deg': 80, 'obsMv': 4.5, 'beta': 0.4},
    {'FOV_deg': 20, 'obsMv': 3.5, 'beta': 0.4},
    {'FOV_deg': 40, 'obsMv': 3.5, 'beta': 0.4},
    {'FOV_deg': 60, 'obsMv': 3.5, 'beta': 0.4},
    {'FOV_deg': 80, 'obsMv': 3.5, 'beta': 0.4},
    {'FOV_deg': 20, 'obsMv': 5.5, 'beta': 0.6},
    {'FOV_deg': 40, 'obsMv': 5.5, 'beta': 0.6},
    {'FOV_deg': 60, 'obsMv': 5.5, 'beta': 0.6},
    {'FOV_deg': 80, 'obsMv': 5.5, 'beta': 0.6},
    {'FOV_deg': 20, 'obsMv': 4.5, 'beta': 0.6},
    {'FOV_deg': 40, 'obsMv': 4.5, 'beta': 0.6},
    {'FOV_deg': 60, 'obsMv': 4.5, 'beta': 0.6},
    {'FOV_deg': 80, 'obsMv': 4.5, 'beta': 0.6},
    {'FOV_deg': 20, 'obsMv': 3.5, 'beta': 0.6},
    {'FOV_deg': 40, 'obsMv': 3.5, 'beta': 0.6},
    {'FOV_deg': 60, 'obsMv': 3.5, 'beta': 0.6},
    {'FOV_deg': 80, 'obsMv': 3.5, 'beta': 0.6},
    {'FOV_deg': 20, 'obsMv': 5.5, 'beta': 0.8},
    {'FOV_deg': 40, 'obsMv': 5.5, 'beta': 0.8},
    {'FOV_deg': 60, 'obsMv': 5.5, 'beta': 0.8},
    {'FOV_deg': 80, 'obsMv': 5.5, 'beta': 0.8},
    {'FOV_deg': 20, 'obsMv': 4.5, 'beta': 0.8},
    {'FOV_deg': 40, 'obsMv': 4.5, 'beta': 0.8},
    {'FOV_deg': 60, 'obsMv': 4.5, 'beta': 0.8},
    {'FOV_deg': 80, 'obsMv': 4.5, 'beta': 0.8},
    {'FOV_deg': 20, 'obsMv': 3.5, 'beta': 0.8},
    {'FOV_deg': 40, 'obsMv': 3.5, 'beta': 0.8},
    {'FOV_deg': 60, 'obsMv': 3.5, 'beta': 0.8},
    {'FOV_deg': 80, 'obsMv': 3.5, 'beta': 0.8},
]

### Parameter ###
# sys
log_dir = './log/'
pattern = patterns[args.i]
# config
seed = 100
roopN = 10000
U = 4096
limitMv = 5.5

FOV_deg = pattern['FOV_deg']
FOV = deg2rad(FOV_deg)
obsMv = pattern['obsMv']
cover_beta = pattern['beta']

sigma = np.arctan2(2*np.tan(FOV*0.5), U)
epsilon = 2*sigma
theta_max = 2*np.arctan2(np.sqrt(2)*np.tan(FOV*0.5), 1)


def main():
    header = [
        'seed', 'N_obs', 'N_candi', 'match', 'time']
    logger = Logger(log_dir, header)
    logger.reset()

    ### Catalog setup ###
    yale = YaleStarCatalog(log_dir=log_dir)
    star_ctlg = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), limitMv)
    star_ctlg.filtering_by_multiple_stars(epsilon)
    #
    angle_ctlg = InterAngleCatalog(
        star_ctlg.get_ID(), star_ctlg.get_RA(), star_ctlg.get_DE(),
        log_dir=log_dir)
    angle_ctlg.create_catalog(theta_max=theta_max, use_log=False)

    ### MONTE CARLO SIMULATION ###
    # observe catalog
    obs_star_ctlg = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), limitMv)
    obs_star_ctlg.filtering_by_multiple_stars(epsilon)
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), obsMv)
    s_vec = equatorial2vec(obs_star_ctlg.get_RA(), obs_star_ctlg.get_DE())
    N_catalog = len(s_vec)

    # preprocess
    tan_FOV2 = np.tan(FOV/2.0)

    # main loop
    for round in tqdm.tqdm(range(roopN)):
        # seed
        sim_seed = gen_seeds(round+seed, 1)[0]
        np.random.seed(seed=sim_seed)
        ### Observed stars ###
        # rotation
        R = special_ortho_group.rvs(dim=3, size=1)
        s_vec_R = np.dot(R, s_vec.T).T
        # in FOV
        cond_1 = s_vec_R[:, 2] > 0.0
        temp = tan_FOV2*s_vec_R[:, 2]
        cond_2 = np.abs(s_vec_R[:, 0]) < temp
        cond_3 = np.abs(s_vec_R[:, 1]) < temp
        in_FOV = cond_1 * cond_2 * cond_3
        # cover rate
        notcover = np.random.binomial(size=N_catalog, n=1, p=1-cover_beta)
        # select stars
        condition = in_FOV*notcover == 1
        obs_ID = obs_star_ctlg.get_ID()[condition]
        obs_RA = obs_star_ctlg.get_RA()[condition]
        obs_DE = obs_star_ctlg.get_DE()[condition]
        ### Add noise ###
        obs_noise = limit_uniform_spherical_vector(
            sigma=sigma, n=len(obs_ID), seed=100)
        axis_vec = np.concatenate(
            [np.cos(obs_RA+0.5*np.pi)[:, np.newaxis],
             np.sin(obs_RA+0.5*np.pi)[:, np.newaxis],
             np.zeros_like(obs_ID)[:, np.newaxis]],
            axis=1)
        R_s_vec = rodrigues_rotation_matrix(0.5*np.pi-obs_DE, axis_vec)
        if len(obs_ID)==1:
            obs_s_vec = rotate_(R_s_vec, obs_noise[np.newaxis, :])
        else:
            obs_s_vec = rotate_(R_s_vec, obs_noise)
        ### Matching ###
        N_obs = len(obs_s_vec)
        if N_obs < 2:
            N_candi = 0
            match = False
        elif N_obs < 5:
            candi_setid_each_list, time_list = matching_set_for_analysis(
                N_obs, obs_s_vec, epsilon,
                star_ctlg.get_RA(), star_ctlg.get_DE(),
                angle_ctlg.get_pair_index(), angle_ctlg.get_inter_angles())
            candi_index = candi_setid_each_list[N_obs - 2]
            N_candi = len(candi_index)
            if N_candi == 1:
                obs_set = set(obs_ID)
                candi_set = set(star_ctlg.get_ID()[candi_index[0]])
                match = obs_set == candi_set
            else:
                match = False
            #
            if N_obs == 2:
                time = time_list[1] - time_list[0]
            else:
                time = time_list[N_obs-1] - time_list[1]
            #
            if N_candi == 0:
                print(f'!!! Warning : N_candi is zero in {round}!!!')
        else:
            match = False
            for di in range(1, N_obs - 3 + 1):
                for dj in range(1, N_obs - di - 2 + 1):
                    for dk in range(1, N_obs - di - dj - 1 + 1):
                        for i in range(N_obs - di - dj - dk):
                            j = i + di
                            k = j + dj
                            l = k + dk
                            candi_setid_each_list, time_list = matching_set_for_analysis(
                                4, obs_s_vec[[i, j, k, l]], epsilon,
                                star_ctlg.get_RA(), star_ctlg.get_DE(),
                                angle_ctlg.get_pair_index(), angle_ctlg.get_inter_angles())
                            candi_index = candi_setid_each_list[4 - 2]
                            N_candi = len(candi_index)
                            if N_candi == 1:
                                obs_set = set(obs_ID[[i, j, k, l]])
                                candi_set = set(star_ctlg.get_ID()[candi_index[0]])
                                match = obs_set == candi_set
                                time = time_list[4-1] - time_list[1]
                            if match:
                                break
                        if match:
                            break
                    if match:
                        break
                if match:
                    break
        ### Log ###
        logger.append([sim_seed, N_obs, N_candi, int(match), time])
    #
    fname = f'stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}'
    logger.save(fname)


def equatorial2vec(alpha, delta):
    rets = np.empty((len(alpha), 3))
    rets[:, 0] = np.cos(alpha) * np.cos(delta)
    rets[:, 1] = np.sin(alpha) * np.cos(delta)
    rets[:, 2] = np.sin(delta)
    return rets


def rotate(R, vec):
    ret = np.empty((R.shape[0], vec.shape[0], vec.shape[1]))
    for i in range(R.shape[0]):
        ret[i, :, :] = np.dot(R[i], vec.T).T
    return ret


def rotate_(R, vec):
    ret = np.empty((vec.shape[0], vec.shape[1]))
    for i in range(R.shape[0]):
        ret[i, :] = np.dot(R[i], vec[i].T).T
    return ret


class Logger:
    def __init__(self, dir, header):
        self.dir = dir
        self.header = header
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        # render
        self.fig = None

    def reset(self):
        self.log = []

    def append(self, l):
        self.log.append(l)

    def save(self, fname):
        df = pd.DataFrame(
            self.log, columns=self.header)
        df.to_csv(self.dir + '/' + f'{fname}.csv')

    def get_df(self):
        df = pd.DataFrame(
            self.log, columns=self.header)
        return df

    def get_array(self):
        array = np.array(self.log)
        return array

    def render(self, FIG_SCALE=3, DPI=100, fps=100):
        if self.fig is None:
            font_setting()
            # fig
            self.N_col = len(self.header)
            ax_col_num = 2
            ax_row_num = math.ceil(self.N_col/2)
            FIG_SIZE = (1.6*FIG_SCALE*ax_col_num, 1*FIG_SCALE*ax_row_num)
            self.fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
            # ax
            self.axes = []
            for i in range(self.N_col):
                ax = self.fig.add_subplot(ax_row_num, ax_col_num, i+1)
                self.axes.append(ax)
        #
        for i in range(self.N_col):
            self.axes[i].clear()
            # self.axes[i].set_xlabel('Count')
            self.axes[i].set_ylabel(self.header[i])
            line = self.axes[i].plot(self.get_array()[:, i], color='blue')
        self.fig.tight_layout()
        #
        plt.pause(0.01)


def analy_result():
    font_setting()
    DPI = 100
    FIG_SIZE = (4, 2.5)
    # 
    logs = []
    header = ['FOV_deg', 'obsMv', 'cover_beta',
             'match_prob', 'match_prob_in_2',
             'match_prob_in_3', 'match_prob_in_4']
    for pattern in patterns:
        FOV_deg = pattern['FOV_deg']
        obsMv = pattern['obsMv']
        cover_beta = pattern['beta']
        # load
        fname = f'stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        # matching prob
        match_prob = df['match'].mean()
        match_prob_in_2 = df[df['N_obs']==2]['match'].mean()
        match_prob_in_3 = df[df['N_obs']==3]['match'].mean()
        match_prob_in_4 = df[df['N_obs']>3]['match'].mean()
        logs.append(
            [FOV_deg, obsMv, cover_beta,
             match_prob, match_prob_in_2,
             match_prob_in_3, match_prob_in_4])
        # obs prob
        fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
        ax = fig.add_subplot(111)
        pm, bins = np.histogram(
            df['N_obs'], bins=range(0, np.max(df['N_obs'])+1), density=True)
        ax.plot(bins[:-1], pm, color='black')
        ax.set_ylabel('Probability mass')
        ax.set_xlabel('The number of stars in FOV')
        ax.set_title(f'Observe probability (FOV:{FOV_deg}, obsMv:{obsMv}, beta:{cover_beta})')
        fig.savefig(f"{log_dir}/inFOV_from_stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.pdf")
    df = pd.DataFrame(logs, columns=header)
    df.to_csv(log_dir+'stats_analy_result.csv')



if __name__ == '__main__':
    main()
    # analy_result()
    print('Task complete')
