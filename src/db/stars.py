import os
import itertools
import numpy as np

from utils.interangle import inter_star_angle_RADE

BASE_DIR = os.path.dirname(__file__)


class StarDB(object):
    def __init__(self, HR, RA, DE, Vmag, Vmax=5.5):
        # filtered by magnitude
        filter = Vmag < Vmax
        self.HR_V = HR[filter]
        self.RA_V = RA[filter]
        self.DE_V = DE[filter]
        self.Vmag_V = Vmag[filter]
        self.N_V = len(self.HR_V)
        self.I_V = np.arange(self.N_V)
        # save
        self.save_cashe()

    def info(self):
        print('----- Star DB -----')
        print(f"the number of stars : {self.N_V}")
        print(
            f"the range of magnitude : [{min(self.Vmag_V)}, {max(self.Vmag_V)}]")

    def get_I(self):
        return self.I_V

    def get_HR(self):
        return self.HR_V.astype('int64')

    def get_RA(self):
        return self.RA_V

    def get_DE(self):
        return self.DE_V

    def get_Vmag(self):
        return self.Vmag_V

    def save_cashe(self):
        np.savetxt(BASE_DIR + '/cache/starDB_I.csv', self.I_V, delimiter=',')
        np.savetxt(BASE_DIR + '/cache/starDB_HR.csv', self.HR_V, delimiter=',')
        np.savetxt(BASE_DIR + '/cache/starDB_RA.csv', self.RA_V, delimiter=',')
        np.savetxt(BASE_DIR + '/cache/starDB_DE.csv', self.DE_V, delimiter=',')
        np.savetxt(
            BASE_DIR + '/cache/starDB_Vmag.csv', self.Vmag_V, delimiter=',')


class StarDBCashe(StarDB):
    def __init__(self):
        self.load_cashe()

    def load_cashe(self):
        self.ID_V = np.loadtxt(
            BASE_DIR + '/cache/starDB_ID.csv', delimiter=',')
        self.RA_V = np.loadtxt(
            BASE_DIR + '/cache/starDB_RA.csv', delimiter=',')
        self.DE_V = np.loadtxt(
            BASE_DIR + '/cache/starDB_DE.csv', delimiter=',')
        self.Vmag_V = np.loadtxt(
            BASE_DIR + '/cache/starDB_Vmag.csv', delimiter=',')
        self.N_V = len(self.ID_V)
        self.I_V = range(self.N_V)
