import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import load_boston


# Hyperparameter
sampling_type = 0
sample_N = 10000
seed = 10
n_obs_list = [2, 3, 4, 5, 6]
Vmax_list = [i+0.5 for i in [1, 2, 3, 4, 5]]
theta_FOV_list = [i*np.pi/180 for i in [5, 10, 30, 60]]
theta_img_list = [np.pi/180*10**i for i in [
    -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]]
k_list = [2.0**i for i in [-1.0, -0.5, 0.0, 0.5, 1.0]]
# system
log_dir = './log/subgraph_monte/'

#
df = pd.read_csv(log_dir + 'summary.csv', index_col=0)
# condition = (df['Vmax'] == 5.5) & (df['theta_FOV']
#                                    == 30*np.pi/180) & (df['k'] == 2**0.5)
# df = df[condition]
# #
df = df.drop('calc_num', axis=1)
# df = df.drop('Vmax', axis=1)
# df = df.drop('theta_FOV', axis=1)
# df = df.drop('k', axis=1)
#
cols = df.columns
rows = cols
corr = df.corr()


def onclick(event):
    x = event.xdata
    y = event.ydata
    axes = event.inaxes
    # axis外クリック時は無視
    if x == None or y == None:
        return
    # ヒートマップ以外をクリック時は無視
    if axes != ax[0]:
        return

    # x、ｙの小数点切り捨ての値がインデックス番号と一致する
    col = cols[int(x)]
    row = rows[int(y)]

    # 散布図の範囲指定用に描くデータの最大値-最小値を計算しておく
    xwid = df[col].max() - df[col].min()
    ywid = df[row].max() - df[row].min()

    # 右の散布図の更新
    ax[1].set_title('Correlation coefficient ' + col + ' - ' + row + ': '
                    + str(round(corr.loc[row, col], 3)))
    ax[1].set_xlabel(col)
    ax[1].set_ylabel(row)
    ln.set_data(df[col], df[row])
    # 散布図の範囲はmax,min値より少し大きめにする
    ax[1].set_xlim(df[col].min() - xwid * 0.1, df[col].max() + xwid * 0.1)
    ax[1].set_ylim(df[row].min() - ywid * 0.1, df[row].max() + ywid * 0.1)
    plt.draw()


fig, ax = plt.subplots(1, 2, figsize=(12, 6))
sns.heatmap(corr, ax=ax[0], vmin=-1, vmax=1,
            annot=True, fmt="3.2f", annot_kws={"size": 5})
ln, = ax[1].plot([np.nan], [np.nan], 'o')  # 散布図プロット
ax[0].set_aspect('equal')
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()
