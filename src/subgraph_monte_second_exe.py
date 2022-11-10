import numpy as np
import time
import subprocess

import psutil


def prepare_args():
    args_list = [
        (0, 2, 0, 0), (0, 2, 0, 1), (0, 3, 0, 0), (0, 3, 0, 1), 
        (1, 0, 0, 0), (1, 0, 0, 1), (1, 1, 0, 0), (1, 1, 0, 1), 
        (1, 2, 0, 0), (1, 2, 0, 1), (1, 3, 0, 0), (1, 3, 0, 1),
        (2, 0, 0, 0), (2, 0, 0, 1), (2, 1, 0, 0), (2, 1, 0, 1), 
        (2, 2, 0, 0), (2, 2, 0, 1), (2, 3, 0, 0), (2, 3, 0, 1), 
        (2, 3, 0, 2), (3, 0, 0, 0), (3, 0, 0, 1), (3, 1, 0, 0), 
        (3, 1, 0, 1), (3, 2, 0, 0), (3, 2, 0, 1), (3, 2, 0, 2), 
        (3, 3, 0, 0), (3, 3, 0, 1), (3, 3, 0, 2), (4, 0, 0, 0), 
        (4, 0, 0, 1), (4, 1, 0, 0), (4, 1, 0, 1), (4, 2, 0, 0), 
        (4, 2, 0, 1), (4, 2, 0, 2), (4, 2, 0, 3), (4, 2, 1, 0), 
        (4, 2, 1, 1), (4, 3, 0, 0), (4, 3, 0, 1), (4, 3, 0, 2), 
        (4, 3, 0, 3), (4, 3, 1, 0), (4, 3, 1, 1), (4, 3, 1, 2), 
        (4, 3, 1, 3), (4, 3, 2, 0), (4, 3, 2, 1), (4, 3, 2, 2),
         (4, 3, 2, 3), (4, 3, 3, 0)
    ]
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
            cmd = f'start python subgraph_monte.py --Vmax {args[0]} --theta_FOV {args[1]} --theta_img {args[2]} --k {args[3]}'
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
