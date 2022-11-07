import numpy as np
import time
import subprocess

import psutil


def prepare_args():
    ### Hyperparameter ###
    # DB
    Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
    theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
    # simulation
    theta_img_list = [np.pi/180*10**i for i in [
        -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
    # subgraph matching
    k_list = [2.0**i for i in [-0.5, 0.0, 0.5, 1.0]]

    ###  ###
    theta_img = len(theta_img_list)
    k = len(k_list)
    Vmax = len(Vmax_list)
    theta_FOV = len(theta_FOV_list)
    #
    args_list = []
    for a in range(theta_img-1, -1, -1):
        for b in range(k-1, -1, -1):
            for c in range(Vmax-1, -1, -1):
                for d in range(theta_FOV-1, -1, -1):
                    args_list.append([a, b, c, d])
    return args_list


def main():
    args_list = prepare_args()
    # init
    counter = 0
    #
    while True:
        cpu = psutil.cpu_percent(interval=1)
        if cpu < 95.0:
            args = args_list[counter]
            cmd = f'start python subgraph_monte.py --theta_img {args[0]} --k {args[1]} --Vmax {args[2]} --theta_FOV {args[3]}'
            res = subprocess.call(cmd, shell=True)
            if res != 0:
                print(args)
            counter += 1
            if counter >= len(args_list):
                break
        time.sleep(1)


if __name__ == '__main__':
    main()
    print('END')
