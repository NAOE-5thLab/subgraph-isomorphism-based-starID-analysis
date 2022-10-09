"""
@author: Kouki Wakita
"""
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils.font import font_setting


class Logger:
    def __init__(self, dir, header):
        self.dir = dir
        self.header = header
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        # render
        self.fig = None

    def reset(self):
        self.log = []

    def append(self, l):
        self.log.append(l)

    def save(self, fname):
        df = pd.DataFrame(
            self.log, columns=self.header)
        df.to_csv(self.dir + '/' + f'{fname}.csv')

    def get_df(self):
        df = pd.DataFrame(
            self.log, columns=self.header)
        return df

    def get_array(self):
        array = np.array(self.log)
        return array

    def render(self, FIG_SCALE=3, DPI=100, fps=100):
        if self.fig is None:
            font_setting()
            # fig
            self.N_col = len(self.header)
            ax_col_num = 2
            ax_row_num = math.ceil(self.N_col/2)
            FIG_SIZE = (1.6*FIG_SCALE*ax_col_num, 1*FIG_SCALE*ax_row_num)
            self.fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
            # ax
            self.axes = []
            for i in range(self.N_col):
                ax = self.fig.add_subplot(ax_row_num, ax_col_num, i+1)
                self.axes.append(ax)
        #
        for i in range(self.N_col):
            self.axes[i].clear()
            # self.axes[i].set_xlabel('Count')
            self.axes[i].set_ylabel(self.header[i])
            line = self.axes[i].plot(self.get_array()[:, i], color='blue')
        self.fig.tight_layout()
        #
        plt.pause(0.01)
