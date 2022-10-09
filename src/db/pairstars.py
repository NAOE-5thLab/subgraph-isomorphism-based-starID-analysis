import os
import itertools
from tkinter import RAISED
import numpy as np

from utils.interangle import inter_star_angle_RADE

BASE_DIR = os.path.dirname(__file__)


class IntervalSearch:
    def __init__(self, y_vector):
        self.y_vector = y_vector
        self.i_vector = np.argsort(self.y_vector)
        self.s_vector = self.y_vector[self.i_vector]

    def get_index_in_certain_interval(self, y_lower_bound, y_upper_bound):
        lower_erase = y_lower_bound < self.s_vector
        upper_erase = self.s_vector < y_upper_bound
        clip = lower_erase*upper_erase
        index = self.i_vector[clip]
        return index


class PairStarDB(object):
    def __init__(self, I, RA, DE, theta_FOV=30*np.pi/180):
        self.RA = RA
        self.DE = DE
        self.create_catalog(I, RA, DE, theta_FOV)
        self.save_cashe()
        #
        self.inr_srch = IntervalSearch(self.Theta_pair_FOV)

    def create_catalog(self, I, RA, DE, theta_FOV):
        ### process for pair stars ###
        I_pair = np.array(list(itertools.combinations(I, 2)))
        Theta_pair = inter_star_angle_RADE(
            RA[I_pair[:, 0]], DE[I_pair[:, 0]], RA[I_pair[:, 1]], DE[I_pair[:, 1]])
        # filtered by FOV
        filter_of_FOV = Theta_pair < 2*theta_FOV
        self.I_pair_FOV = I_pair[filter_of_FOV]
        self.Theta_pair_FOV = Theta_pair[filter_of_FOV]

    def info(self):
        print('----- Pair Star DB -----')
        print(f"the number of pairs : {len(self.I_pair_FOV)}")
        print(
            f"the range of inter angle [deg.] : [{min(self.Theta_pair_FOV)*180/np.pi}, {max(self.Theta_pair_FOV)*180/np.pi}]")

    def get_RA(self):
        return self.RA

    def get_DE(self):
        return self.DE

    def get_I_pair_FOV(self):
        return self.I_pair_FOV

    def get_Theta_pair_FOV(self):
        return self.Theta_pair_FOV

    def get_pairIDs_in_certain_interval(self, y_lower_bound, y_upper_bound):
        index = self.inr_srch.get_index_in_certain_interval(
            y_lower_bound, y_upper_bound)
        return self.I_pair_FOV[index]

    def get_index_in_certain_interval(self, y_lower_bound, y_upper_bound):
        return self.inr_srch.get_index_in_certain_interval(y_lower_bound, y_upper_bound)

    def save_cashe(self):
        np.savetxt(BASE_DIR + '/cache/I_pair_FOV.csv',
                   self.I_pair_FOV, delimiter=',')
        np.savetxt(BASE_DIR + '/cache/Theta_pair_FOV.csv',
                   self.Theta_pair_FOV, delimiter=',')


class PairStarDBCashe(PairStarDB):
    def __init__(self):
        self.load_cashe()
        self.inr_srch = IntervalSearch(self.Theta_pair_FOV)

    def load_cashe(self):
        self.I_pair_FOV = np.loadtxt(
            BASE_DIR + '/cache/I_pair_FOV.csv', delimiter=',')
        self.Theta_pair_FOV = np.loadtxt(
            BASE_DIR + '/cache/Theta_pair_FOV.csv', delimiter=',')
