import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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


def main():
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
    # plot
    font_setting()
    DPI = 200
    FIG_SIZE = (8, 3.5)
    marker = 'o'

    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=4)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('$P(N_{\mathrm{obs}} >= 4 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'Nobs_le_4.pdf'
    fig.savefig(log_dir+fname)


    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['P(n>=10)'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('$P(N_{\mathrm{obs}} >= 10 \mid M_v, \\theta_{\mathrm{FOV}}, \\beta)$')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'Nobs_le_10.pdf'
    fig.savefig(log_dir+fname)

    # 
    fname = 'stats_analy_result.csv'
    df = pd.read_csv(log_dir+fname)

    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    df_temp = df[df['obsMv'] == 3.5][df['cover_beta'] == 0.0]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 3.5][df['cover_beta'] == 0.2]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 3.5][df['cover_beta'] == 0.4]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 3.5][df['cover_beta'] == 0.6]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 3.5][df['cover_beta'] == 0.8]
    ax1.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax1.set_title('$M_{v} = 3.5$')
    ax1.set_xlabel('FOV [deg.]')
    ax1.set_ylabel('Matching Probability')
    ax1.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 4.5][df['cover_beta'] == 0.0]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 4.5][df['cover_beta'] == 0.2]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 4.5][df['cover_beta'] == 0.4]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 4.5][df['cover_beta'] == 0.6]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 4.5][df['cover_beta'] == 0.8]
    ax2.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax2.set_title('$M_{v} = 4.5$')
    ax2.set_xlabel('FOV [deg.]')
    # ax2.set_ylabel('Matching Probability')
    ax2.set_ylim(-0.05, 1.05)

    df_temp = df[df['obsMv'] == 5.5][df['cover_beta'] == 0.0]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(0/5), label='$\\beta = 0.0$')
    df_temp = df[df['obsMv'] == 5.5][df['cover_beta'] == 0.2]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(1/5), label='$\\beta = 0.2$')
    df_temp = df[df['obsMv'] == 5.5][df['cover_beta'] == 0.4]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(2/5), label='$\\beta = 0.4$')
    df_temp = df[df['obsMv'] == 5.5][df['cover_beta'] == 0.6]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(3/5), label='$\\beta = 0.6$')
    df_temp = df[df['obsMv'] == 5.5][df['cover_beta'] == 0.8]
    ax3.plot(df_temp['FOV_deg'], df_temp['match_prob'], marker=marker, color = cm.viridis(4/5), label='$\\beta = 0.8$')
    ax3.set_title('$M_{v} = 5.5$')
    ax3.set_xlabel('FOV [deg.]')
    # ax3.set_ylabel('Matching Probability')
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend()

    fig.tight_layout()
    fname = 'matching_prob.pdf'
    fig.savefig(log_dir+fname)


if __name__ == '__main__':
    main()
    print('Task complete')

