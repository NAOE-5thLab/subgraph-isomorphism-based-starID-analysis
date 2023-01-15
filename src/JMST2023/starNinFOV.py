import os
from line_profiler import LineProfiler
import argparse

import tqdm
import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import special_ortho_group
from scipy.stats import binom

from database import YaleStarCatalog, StarCatalog
from utils.unit_convert import deg2rad
from utils.font import font_setting

parser = argparse.ArgumentParser()
parser.add_argument(
    "--i",
    type=int,
    default=0,
    help="degree")
args = parser.parse_args()

patterns = [
    {'FOV_deg': 20, 'obsMv': 5.5, 'alpha': 0.0},
    {'FOV_deg': 40, 'obsMv': 5.5, 'alpha': 0.0},
    {'FOV_deg': 60, 'obsMv': 5.5, 'alpha': 0.0},
    {'FOV_deg': 80, 'obsMv': 5.5, 'alpha': 0.0},
    {'FOV_deg': 20, 'obsMv': 4.5, 'alpha': 0.0},
    {'FOV_deg': 40, 'obsMv': 4.5, 'alpha': 0.0},
    {'FOV_deg': 60, 'obsMv': 4.5, 'alpha': 0.0},
    {'FOV_deg': 80, 'obsMv': 4.5, 'alpha': 0.0},
    {'FOV_deg': 20, 'obsMv': 3.5, 'alpha': 0.0},
    {'FOV_deg': 40, 'obsMv': 3.5, 'alpha': 0.0},
    {'FOV_deg': 60, 'obsMv': 3.5, 'alpha': 0.0},
    {'FOV_deg': 80, 'obsMv': 3.5, 'alpha': 0.0},
    {'FOV_deg': 20, 'obsMv': 5.5, 'alpha': 0.2},
    {'FOV_deg': 40, 'obsMv': 5.5, 'alpha': 0.2},
    {'FOV_deg': 60, 'obsMv': 5.5, 'alpha': 0.2},
    {'FOV_deg': 80, 'obsMv': 5.5, 'alpha': 0.2},
    {'FOV_deg': 20, 'obsMv': 4.5, 'alpha': 0.2},
    {'FOV_deg': 40, 'obsMv': 4.5, 'alpha': 0.2},
    {'FOV_deg': 60, 'obsMv': 4.5, 'alpha': 0.2},
    {'FOV_deg': 80, 'obsMv': 4.5, 'alpha': 0.2},
    {'FOV_deg': 20, 'obsMv': 3.5, 'alpha': 0.2},
    {'FOV_deg': 40, 'obsMv': 3.5, 'alpha': 0.2},
    {'FOV_deg': 60, 'obsMv': 3.5, 'alpha': 0.2},
    {'FOV_deg': 80, 'obsMv': 3.5, 'alpha': 0.2},
    {'FOV_deg': 20, 'obsMv': 5.5, 'alpha': 0.4},
    {'FOV_deg': 40, 'obsMv': 5.5, 'alpha': 0.4},
    {'FOV_deg': 60, 'obsMv': 5.5, 'alpha': 0.4},
    {'FOV_deg': 80, 'obsMv': 5.5, 'alpha': 0.4},
    {'FOV_deg': 20, 'obsMv': 4.5, 'alpha': 0.4},
    {'FOV_deg': 40, 'obsMv': 4.5, 'alpha': 0.4},
    {'FOV_deg': 60, 'obsMv': 4.5, 'alpha': 0.4},
    {'FOV_deg': 80, 'obsMv': 4.5, 'alpha': 0.4},
    {'FOV_deg': 20, 'obsMv': 3.5, 'alpha': 0.4},
    {'FOV_deg': 40, 'obsMv': 3.5, 'alpha': 0.4},
    {'FOV_deg': 60, 'obsMv': 3.5, 'alpha': 0.4},
    {'FOV_deg': 80, 'obsMv': 3.5, 'alpha': 0.4},
    {'FOV_deg': 20, 'obsMv': 5.5, 'alpha': 0.6},
    {'FOV_deg': 40, 'obsMv': 5.5, 'alpha': 0.6},
    {'FOV_deg': 60, 'obsMv': 5.5, 'alpha': 0.6},
    {'FOV_deg': 80, 'obsMv': 5.5, 'alpha': 0.6},
    {'FOV_deg': 20, 'obsMv': 4.5, 'alpha': 0.6},
    {'FOV_deg': 40, 'obsMv': 4.5, 'alpha': 0.6},
    {'FOV_deg': 60, 'obsMv': 4.5, 'alpha': 0.6},
    {'FOV_deg': 80, 'obsMv': 4.5, 'alpha': 0.6},
    {'FOV_deg': 20, 'obsMv': 3.5, 'alpha': 0.6},
    {'FOV_deg': 40, 'obsMv': 3.5, 'alpha': 0.6},
    {'FOV_deg': 60, 'obsMv': 3.5, 'alpha': 0.6},
    {'FOV_deg': 80, 'obsMv': 3.5, 'alpha': 0.6},
    {'FOV_deg': 20, 'obsMv': 5.5, 'alpha': 0.8},
    {'FOV_deg': 40, 'obsMv': 5.5, 'alpha': 0.8},
    {'FOV_deg': 60, 'obsMv': 5.5, 'alpha': 0.8},
    {'FOV_deg': 80, 'obsMv': 5.5, 'alpha': 0.8},
    {'FOV_deg': 20, 'obsMv': 4.5, 'alpha': 0.8},
    {'FOV_deg': 40, 'obsMv': 4.5, 'alpha': 0.8},
    {'FOV_deg': 60, 'obsMv': 4.5, 'alpha': 0.8},
    {'FOV_deg': 80, 'obsMv': 4.5, 'alpha': 0.8},
    {'FOV_deg': 20, 'obsMv': 3.5, 'alpha': 0.8},
    {'FOV_deg': 40, 'obsMv': 3.5, 'alpha': 0.8},
    {'FOV_deg': 60, 'obsMv': 3.5, 'alpha': 0.8},
    {'FOV_deg': 80, 'obsMv': 3.5, 'alpha': 0.8},
]

### Parameter ###
# sys
log_dir = './log/'
pattern = patterns[args.i]
# config
seed = 100
roopN = int(1e5)
calcN = int(1e4)
U = 4096
limitMv = 5.5

FOV_deg = pattern['FOV_deg']
FOV = deg2rad(FOV_deg)
obsMv = pattern['obsMv']
cover_alpha = pattern['alpha']

sigma = np.arctan2(2*np.tan(FOV*0.5), U)
epsilon = 2*sigma
theta_max = 2*np.arctan2(np.sqrt(2)*np.tan(FOV*0.5), 1)


def simulation():
    # check file
    path = f'{log_dir}/inFOV_{FOV_deg}_{obsMv}.dat'
    if os.path.isfile(path):
        print(f'{path} exist!')
        return 0
    
    ### Catalog setup ###
    yale = YaleStarCatalog(log_dir=log_dir)
    # 
    obs_star_ctlg = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), limitMv)
    obs_star_ctlg.filtering_by_multiple_stars(epsilon)
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), obsMv)
    s_vec = equatorial2vec(obs_star_ctlg.get_RA(), obs_star_ctlg.get_DE())

    ### Monte Carlo ###
    tan_FOV2 = np.tan(FOV/2.0)
    # loop
    log_file = open(path, 'w')
    for round in tqdm.tqdm(range(roopN)):
        # sampling
        np.random.seed(seed=round+seed)
        R = special_ortho_group.rvs(dim=3, size=calcN)
        eq_vec_R = rotate(R, s_vec)
        # FOV condition
        cond_1 = eq_vec_R[:, :, 2] > 0.0
        temp = tan_FOV2*eq_vec_R[:, :, 2]
        cond_2 = np.abs(eq_vec_R[:, :, 0]) < temp
        cond_3 = np.abs(eq_vec_R[:, :, 1]) < temp
        in_FOV = cond_1 * cond_2 * cond_3
        in_FOV_n = np.sum(in_FOV, axis=1)
        # save
        for j in range(calcN):
            log_file.write(f'{in_FOV_n[j]}\n')
    log_file.close()
    return 0


@jit(nopython=True)
def equatorial2vec(alpha, delta):
    rets = np.empty((len(alpha), 3))
    rets[:, 0] = np.cos(alpha) * np.cos(delta)
    rets[:, 1] = np.sin(alpha) * np.cos(delta)
    rets[:, 2] = np.sin(delta)
    return rets


@jit(nopython=True)
def rotate(R, vec):
    ret = np.empty((R.shape[0], vec.shape[0], vec.shape[1]))
    for i in range(R.shape[0]):
        ret[i, :, :] = np.dot(R[i], vec.T).T
    return ret


def analy_result():
    font_setting()
    DPI = 100
    FIG_SIZE = (4, 2.5)
    # 
    max_n = 1000
    true_max_n = 0
    counter = 0
    n_counter = np.zeros(max_n)
    path = f'{log_dir}/inFOV_{FOV_deg}_{obsMv}.dat'
    for line in line_gene(path):
        n = int(line)
        counter += 1
        n_counter[n] += 1
        if true_max_n < n:
            true_max_n = n
    # 
    pNs = n_counter/counter
    pns = np.zeros(max_n)
    for n in range(max_n):
        for N in range(max_n):
            if pNs[N] > 0.0:
                pn_N = binom.pmf(n, N, 1-cover_alpha)
                pns[n] += pNs[N] * pn_N
            else:
                pns[n] += 0.0
    print(pns.sum(), pNs.sum())
    # save
    df = pd.DataFrame(range(true_max_n+1), pNs[:true_max_n+1], columns=['n', 'prob'])
    df.to_csv(f'{log_dir}/inFOV_prob_FOV{FOV_deg}_obsMv{obsMv}.csv')
    df = pd.DataFrame(range(true_max_n+1), pns[:true_max_n+1], columns=['n', 'prob'])
    df.to_csv(f'{log_dir}/inFOV_prob_FOV{FOV_deg}_obsMv{obsMv}_alpha{cover_alpha}.csv')
    # plot 
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax = fig.add_subplot(111)
    ax.plot(
        range(true_max_n+1), pNs[:true_max_n+1],
        color='black',
        label='$P(N^{\prime}_{\mathrm{obs}} \mid M_v, \\theta_{\mathrm{FOV}})$')
    ax.plot(
        range(true_max_n+1), pns[:true_max_n+1],
        color='black',
        label='$P(N_{\mathrm{obs}} \mid M_v, \\theta_{\mathrm{FOV}}), \\beta$')
    ax.set_ylabel('Probability mass')
    ax.set_xlabel('The number of stars in FOV')
    fig.savefig(f"{log_dir}/inFOV_FOV{FOV_deg}_obsMv{obsMv}_alpha{cover_alpha}.pdf")


def line_gene(path):
        with open(path, 'r') as f:
            for line in f:
                yield line
                 

def plot_starnum_prob(star_num):
    font_setting()
    DPI = 100
    FIG_SIZE = (5, 3.5)
    ls_list = ["solid", "dashed", "dashdot"]
    color_list = ['black', 'b', 'r', 'g']
    span = 0.05
    alphas = np.arange(0, 1+span, span)
    # 
    df = pd.DataFrame(alphas, columns=['alpha'])
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            print(f'Magnitude Limit = {Mv_obs} [Mv], FOV = {FOV_deg} [deg.]')
            with open(f"{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_dens.dat", mode='r') as f:
                lines = f.readlines()
            pNs = np.zeros(len(lines))
            count_all = 0
            for line in lines:
                N, count = line.split(',')
                pNs[int(N)] = int(count[:-3])
                count_all += int(count[:-3])
            pNs = pNs/count_all
            # 
            probs = []
            for alpha in alphas:
                pns = np.zeros(len(lines))
                for n in range(len(lines)):
                    for N in range(len(lines)):
                        pn_N = binom.pmf(n, N, 1-alpha)
                        pns[n] += pNs[N] * pn_N
                # 
                prob = 0
                for n in range(star_num+1):
                    prob += pns[n]
                probs.append(prob)
            # 
            df[f'Mv{Mv_obs}_FOV{FOV_deg}'] = probs
    df.to_csv(f"{log_dir}/inFOV_starnum_prob_{star_num}.csv")
    # 
    df = pd.read_csv(f"{log_dir}/inFOV_starnum_prob_{star_num}.csv")
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax = fig.add_subplot(111)
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            probs = df[f'Mv{Mv_obs}_FOV{FOV_deg}']
            ax.plot(
                alphas, probs,
                color=color_list[j],
                linestyle=ls_list[i],
                label=f'{Mv_obs}, {FOV_deg}')
    ax.set_ylabel(f'Probability $P (n \leq {star_num})$')
    ax.set_xlabel('Cover rate $\\alpha$')
    ax.legend(title = "Mag, FOV", bbox_to_anchor=(1.05, 1), loc='upper left')
    fig.tight_layout()
    fig.savefig(f"{log_dir}/inFOV_starnum_prob_{star_num}.pdf")
    
    
def main():
    simulation()
    analy_result()
    # plot_prob(n=4)
    # plot_prob(n=10)

if __name__ == '__main__':
    # prof = LineProfiler()
    # prof.add_function(main)
    # prof.runcall(main)
    # prof.print_stats()
    main()
    print('task complete')
