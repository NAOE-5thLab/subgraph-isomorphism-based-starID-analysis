from re import sub
import time
import subprocess

import psutil


def prepare_args():
    theta_img = 10
    k = 5
    Vmax = 5
    theta_FOV = 4
    #
    args_list = []
    for a in range(theta_img):
        for b in range(k):
            for c in range(Vmax):
                for d in range(theta_FOV):
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

if __name__=='__main__':
    main()
    print('END')
