import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binom

from utils.font import font_setting


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


def analy_result():
    font_setting()
    DPI = 100
    FIG_SIZE = (8, 5)
    #
    for pattern in patterns:
        # param
        FOV_deg = pattern['FOV_deg']
        obsMv = pattern['obsMv']
        cover_beta = pattern['beta']
        # fig
        fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
        ax = fig.add_subplot(111)

        ### Result of starid.py ###
        fname = f'stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        pm, bins = np.histogram(
            df['N_obs'], bins=range(0, np.max(df['N_obs'])+1), density=True)
        ax.plot(
            bins[:-1], pm,
            color='black',
            label='$\hat{P}(N_{\mathrm{obs}} \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')

        ### Result of starid.py ###
        fname = f'inFOV_prob_FOV{FOV_deg}_obsMv{obsMv}.csv'
        df = pd.read_csv(log_dir+fname)
        pNs = df['prob'].to_numpy()
        ax.plot(
            range(len(pNs)), pNs,
            color='red',
            ls='dashed',
            label='$P(N^{\prime}_{\mathrm{obs}} \mid M_v, \\theta_{\mathrm{FOV}})$')
        #
        fname = f'inFOV_prob_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        pns = df['prob'].to_numpy()
        ax.plot(
            range(len(pns)), pns,
            color='red',
            label='$P(N_{\mathrm{obs}} \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')

        ### plot ###
        ax.set_ylabel('Probability mass')
        ax.set_xlabel('The number of stars in FOV')
        ax.set_title(f'FOV:{FOV_deg}, obsMv:{obsMv}, beta:{cover_beta}')
        ax.legend()
        fig.savefig(f"{log_dir}/inFOV_compare_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.pdf")


if __name__ == '__main__':
    # main()
    analy_result()
    print('Task complete')



































