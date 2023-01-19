import os
import time
import subprocess

import psutil


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


def main():
    # init
    counter = 0
    #
    while True:
        cpu = psutil.cpu_percent(interval=1)
        if cpu < 95.0:
            print(f'start {counter}/{len(patterns)} (cpu {cpu}%)')
            cmd = f'start /min python starid.py --i {counter}'
            res = subprocess.call(cmd, shell=True)
            if res != 0:
                print(patterns[counter])
            counter += 1
            if counter >= len(patterns):
                break
        time.sleep(10)


if __name__ == '__main__':
    main()
    print('END')
