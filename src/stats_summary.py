import itertools
import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import Param
import utils


def main():
    utils.font_setting()
    #
    conf = Param()
    ps = list(itertools.product(
        conf.n_obs_list, conf.Vmax_list,
        conf.theta_FOV_list, conf.theta_img_list, conf.k_list))
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
        fname = f'stats_{n_obs}_{Vmax}_{theta_FOV*180/np.pi}_{theta_img*180/np.pi}_{k}'
        df = pd.read_csv(conf.log_dir + fname + '.csv', index_col=0)

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
    df.to_csv(conf.log_dir + f'summary.csv')


if __name__ == '__main__':
    main()
    print('task completed')
