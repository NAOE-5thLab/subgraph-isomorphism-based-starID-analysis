import os
import time
import subprocess

import numpy as np
import psutil

from config import Param


def main():
    conf = Param()
    #
    args_list = []
    # for i_Vmax in range(len(conf.Vmax_list)):
    #     for i_theta_FOV in range(len(conf.theta_FOV_list)):
    #         for i_theta_img in range(len(conf.theta_img_list)):
    #             for i_k in range(len(conf.k_list)):
    #                 for i_parallel in range(conf.parallel_num):
    #                     args_list.append(
    #                         [i_Vmax, i_theta_FOV, i_theta_img, i_k, i_parallel])
    for i_Vmax in range(len(conf.Vmax_list)-1, -1, -1):
        for i_theta_FOV in range(len(conf.theta_FOV_list)-1, -1, -1):
            for i_theta_img in range(len(conf.theta_img_list)-1, -1, -1):
                for i_k in range(len(conf.k_list)-1, -1, -1):
                    for i_parallel in range(conf.parallel_num):
                        args_list.append(
                            [i_Vmax, i_theta_FOV, i_theta_img, i_k, i_parallel])
    # init
    counter = 0
    #
    while True:
        args = args_list[counter]

        fname = f'stats_{conf.n_obs_list[-1]}_{conf.Vmax_list[args[0]]}_{conf.theta_FOV_list[args[1]]*180/np.pi}_{conf.theta_img_list[args[2]]*180/np.pi}_{conf.k_list[args[3]]}'
        if os.path.isfile(conf.log_dir + '/' + f'{fname}.csv'):
            counter += 1
            if counter >= len(args_list):
                break
            else:
                continue
        #
        cpu = psutil.cpu_percent(interval=1)
        print(f'cpu : {cpu}')
        if cpu < 95.0:
            cmd = f'start /min python subgraph_monte.py --i_Vmax {args[0]} --i_theta_FOV {args[1]} --i_theta_img {args[2]} --i_k {args[3]} --i_parallel {args[4]}'
            print(f'start : {args}')
            res = subprocess.call(cmd, shell=True)
            if res != 0:
                print(args)
            counter += 1
            if counter >= len(args_list):
                break
        time.sleep(5)


if __name__ == '__main__':
    main()
    print('END')
