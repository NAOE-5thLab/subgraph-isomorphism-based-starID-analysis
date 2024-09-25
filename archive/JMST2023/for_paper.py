import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from database import YaleStarCatalog, StarCatalog, InterAngleCatalog
from utils.unit_convert import deg2rad
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

DPI = 200
FIG_SIZE = (6, 2.7)
marker = 'o'
font_setting()


def main():
    plot_inFOV()
    plot_matching()
    plot_db()


def plot_inFOV():
    logs = []
    for pattern in patterns:
        # param
        FOV_deg = pattern['FOV_deg']
        obsMv = pattern['obsMv']
        cover_beta = pattern['beta']
        log = pattern
        #
        fname = f'inFOV_prob_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        pns = df['prob'].to_numpy()
        #
        log['P(n>=4)'] = pns[4:].sum()
        log['P(n>=10)'] = pns[10:].sum()
        #
        logs.append(log)
    df = pd.DataFrame(logs)
    ### plot ###
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'Nobs_le_4.pdf'
    fig.savefig(log_dir+fname)

    ### plot ###
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'Nobs_le_10.pdf'
    fig.savefig(log_dir+fname)


def plot_matching():
    logs = []
    for pattern in patterns:
        # param
        FOV_deg = pattern['FOV_deg']
        obsMv = pattern['obsMv']
        cover_beta = pattern['beta']
        log = pattern
        #
        fname = f'stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        match = df['match'].to_numpy()
        #
        log['match_prob'] = match.mean()
        logs.append(log)
    df = pd.DataFrame(logs)

    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('Matching Probability')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('Matching Probability')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, ls='dotted', color=cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('Matching Probability')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'matching_prob.pdf'
    fig.savefig(log_dir+fname)


def plot_db():
    ### db size ###
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    #
    data = []
    for i, FOV_deg in enumerate([20, 40, 60, 80]):
        U = 4096
        limitMv = 5.5
        FOV = deg2rad(FOV_deg)
        sigma = np.arctan2(2*np.tan(FOV*0.5), U)
        epsilon = 2*sigma
        theta_max = 2*np.arctan2(np.sqrt(2)*np.tan(FOV*0.5), 1)
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
        angles = angle_ctlg.get_inter_angles()
        #
        data.append([FOV_deg, len(angles)])
    data = np.array(data)
    ax1.plot(data[:, 0], data[:, 1], marker=marker, ls='dotted', color = cm.viridis(9/10))
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('Database size ($\mathcal{I}_{\mathrm{pair}}$)')

    ### matching time ###
    t_2_data = {20: [], 40: [], 60: [], 80: []}
    t_3_data = {20: [], 40: [], 60: [], 80: []}
    t_4_data = {20: [], 40: [], 60: [], 80: []}
    for pattern in patterns:
        # param
        FOV_deg = pattern['FOV_deg']
        obsMv = pattern['obsMv']
        cover_beta = pattern['beta']
        #
        fname = f'stats_FOV{FOV_deg}_obsMv{obsMv}_beta{cover_beta}.csv'
        df = pd.read_csv(log_dir+fname)
        N_obs = df['N_obs'].to_numpy()
        time = df['time'].to_numpy()
        #
        t_2_data[FOV_deg].append(time[N_obs == 2])
        t_3_data[FOV_deg].append(time[N_obs == 3])
        t_4_data[FOV_deg].append(time[N_obs >= 4])

    logs = []
    for FOV_deg in [20, 40, 60, 80]:
        log = {'FOV_deg': FOV_deg}
        log['t_2'] = np.hstack(t_2_data[FOV_deg]).mean()
        log['t_3'] = np.hstack(t_3_data[FOV_deg]).mean()
        log['t_4'] = np.hstack(t_4_data[FOV_deg]).mean()
        logs.append(log)
    df = pd.DataFrame(logs)

    ax2.plot(df['FOV_deg'], df['t_2'], marker=marker, ls='dotted', color=cm.viridis(0/3), label='$N_{\mathrm{obs}} = 2$')
    ax2.plot(df['FOV_deg'], df['t_3'], marker=marker, ls='dotted', color=cm.viridis(1/3), label='$N_{\mathrm{obs}} = 3$')
    ax2.plot(df['FOV_deg'], df['t_4'], marker=marker, ls='dotted', color=cm.viridis(2/3), label='$N_{\mathrm{obs}} \geq 4$')
    # ax2.set_yscale('log')
    ax2.set_xlabel('FOV [deg.]')
    ax2.set_ylabel('Runtime [s]')
    # ax2.set_ylim(-0.0005, 10.0)
    ax2.legend()

    fig.tight_layout()
    fname = 'runtime.pdf'
    fig.savefig(log_dir+fname)


def convert_for_step_plot(bins, pd):
    x = []
    y = []
    for i in range(len(bins)-1):
        x.append(bins[i])
        y.append(pd[i])
        x.append(bins[i+1])
        y.append(pd[i])
    return x, y

if __name__ == '__main__':
    main()
    print('Task complete')

