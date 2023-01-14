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
    "--FOV",
    type=int,
    default=20,
    help="degree")
parser.add_argument(
    "--Mv",
    type=float,
    default=5.5,
    help="visual magnitude")
args = parser.parse_args()

### Parameter ###
# sys
log_dir = './log/'
# config
master_seed = 10
FOV_deg = args.FOV
FOV = deg2rad(FOV_deg)
roopN = int(1e5)
calcN = int(1e4)
# catalog
Mv_max = 5.5
Mv_obs = args.Mv
multi_angle = 1.0e-5


def main():
    sim()
    # hist_csv()
    # plot_prob_alpha(0.2)
    # plot_prob_alpha(0.4)
    # plot_prob_alpha(0.6)
    # plot_prob_alpha(0.8)
    # plot_starnum_prob(10)
    # plot_all()
    # prob_n(n=4)
    # prob_n(n=10)
    # plot_prob(n=4)
    # plot_prob(n=10)


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


def sim():
    ### Catalog setup ###
    yale = YaleStarCatalog(log_dir=log_dir)
    # 
    obs_star_ctlg = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), Mv_max)
    obs_star_ctlg.filtering_by_multiple_stars(multi_angle)    
    obs_star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), Mv_obs)
    eq_vec = equatorial2vec(obs_star_ctlg.get_RA(), obs_star_ctlg.get_DE())

    ### Monte Carlo ###
    tan_FOV2 = np.tan(FOV/2.0)
    # loop
    log_file = open(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}.dat', 'w')
    for i in tqdm.tqdm(range(roopN)):
        # sampling
        np.random.seed(seed=i+master_seed)
        R = special_ortho_group.rvs(dim=3, size=calcN)
        eq_vec_R = rotate(R, eq_vec)
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


def hist_csv():
    def line_gene(path):
        with open(path, 'r') as f:
            for line in f:
                yield line
    width_list = [60, 180, 630]
    #
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            print(f'FOV : {FOV_deg}, Mv : {Mv_obs}')
            progress = 0
            counter = np.zeros(width_list[i]+1)
            for line in line_gene(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}.dat'):
                num = int(line)
                counter[num] += 1
                progress += 1
                if progress % 1e8==0:
                    print(f'progress : {progress}')
            with open(f"{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_dens.dat", mode='w') as f:
                for k in range(width_list[i]+1):
                    f.write(f'{k},{counter[k]}\n')


def plot_all():
    font_setting()
    DPI = 100
    FIG_SIZE = (4, 2.5)
    #
    width_list = [60, 180, 630]
    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    color_list = ['black', 'b', 'r', 'g']
    #
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
        ax = fig.add_subplot(111)
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            print(f'FOV : {FOV_deg}, Mv : {Mv_obs}')
            print(f'loading data ...')
            log_file = open(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}.dat', 'r')
            num_list = log_file.readlines()
            num_list = [int(num) for num in num_list]
            log_file.close()
            #
            print(f'calc hist ...')
            pd, bins = np.histogram(num_list, bins=range(
                0, width_list[i]+1), density=True)
            print(f'ploting ...')
            ax.plot(
                bins[:-1], pd,
                color=color_list[j],
                linestyle=ls_list[j],
                label=f'FOV = {FOV_deg} [deg.]',
            )
        ax.set_ylabel('Probability')
        ax.set_xlabel('The number of stars in FOV')
        ax.legend()
        fig.savefig(f"{log_dir}/inFOV_Mv{Mv_obs}_dens.pdf")


def prob_n(n=10):
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            print(f'FOV : {FOV_deg}, Mv : {Mv_obs}')
            print(f'loading data ...')
            load_file = open(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}.dat', 'r')
            num_list = load_file.readlines()
            nums = np.array([int(num) for num in num_list])
            load_file.close()
            #
            log_file = open(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_prob_{n}.dat', 'w')
            for alpha in range(0, 1+0.05, 0.05):
                freq = np.sum(nums < n/alpha)
                prob = freq/len(nums)
                log_file.write(f'{alpha},{prob}\n')
            log_file.close()


def plot_prob(n=10):
    font_setting()
    DPI = 100
    FIG_SIZE = (4, 2.5)
    # 
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax = fig.add_subplot(111)
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
            print(f'FOV : {FOV_deg}, Mv : {Mv_obs}')
            print(f'loading data ...')
            load_file = open(f'{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_prob_{n}.dat', 'r')
            alphas = []
            probs = []
            for line in load_file.readlines():
                alphas.append(float(line.split(',')[0]))
                probs.append(float(line.split(',')[1]))
            load_file.close()
            ax.plot(alphas, probs)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Alpha')
    ax.set_xlabel('Probability')
    
    
def plot_prob_alpha(alpha):
    font_setting()
    DPI = 100
    FIG_SIZE = (4, 2.5)
    width_list = [60, 180, 630]
    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    color_list = ['black', 'b', 'r', 'g']
    for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
        fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
        ax = fig.add_subplot(111)
        for j, FOV_deg in enumerate([20, 40, 60, 80]):
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
            pns = np.zeros(len(lines))
            for n in range(len(lines)):
                for N in range(len(lines)):
                    # alpha_prob = alpha**(N-n)*(1-alpha)**n
                    # pn_N = comb(N, n, exact=True)*alpha_prob
                    pn_N = binom.pmf(n, N, 1-alpha)
                    pns[n] += pNs[N] * pn_N
            # 
            print(pns.sum(), pNs.sum())
            #
            ax.plot(
                range(len(pns)), pns,
                color=color_list[j],
                linestyle=ls_list[j],
                label=f'FOV = {FOV_deg} [deg.]',
            )
        # 
        cum_p = 0
        for k, pn in enumerate(pns):
            cum_p += pn
            if 1.0 - cum_p < 1.0e-8:
                if pn == 0.0:
                    print(f'end : {k}')
                    ax.set_xlim(0, k+1)
                    break
        # 
        ax.set_ylabel('Probability')
        ax.set_xlabel('The number of stars in FOV')
        ax.legend()
        fig.savefig(f"{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_{alpha}_dens.pdf")
        

def plot_starnum_prob(star_num):
    font_setting()
    DPI = 100
    FIG_SIZE = (5, 3.5)
    ls_list = ["solid", "dashed", "dashdot"]
    color_list = ['black', 'b', 'r', 'g']
    span = 0.05
    alphas = np.arange(0, 1+span, span)
    # 
    # df = pd.DataFrame(alphas, columns=['alpha'])
    # for i, Mv_obs in enumerate([3.5, 4.5, 5.5]):
    #     for j, FOV_deg in enumerate([20, 40, 60, 80]):
    #         print(f'Magnitude Limit = {Mv_obs} [Mv], FOV = {FOV_deg} [deg.]')
    #         with open(f"{log_dir}/inFOV_{FOV_deg}_{Mv_obs}_dens.dat", mode='r') as f:
    #             lines = f.readlines()
    #         pNs = np.zeros(len(lines))
    #         count_all = 0
    #         for line in lines:
    #             N, count = line.split(',')
    #             pNs[int(N)] = int(count[:-3])
    #             count_all += int(count[:-3])
    #         pNs = pNs/count_all
    #         # 
    #         probs = []
    #         for alpha in alphas:
    #             pns = np.zeros(len(lines))
    #             for n in range(len(lines)):
    #                 for N in range(len(lines)):
    #                     pn_N = binom.pmf(n, N, 1-alpha)
    #                     pns[n] += pNs[N] * pn_N
    #             # 
    #             prob = 0
    #             for n in range(star_num+1):
    #                 prob += pns[n]
    #             probs.append(prob)
    #         # 
    #         df[f'Mv{Mv_obs}_FOV{FOV_deg}'] = probs
    # df.to_csv(f"{log_dir}/inFOV_starnum_prob_{star_num}.csv")
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
    

if __name__ == '__main__':
    # prof = LineProfiler()
    # prof.add_function(main)
    # prof.runcall(main)
    # prof.print_stats()
    main()
    print('task complete')
