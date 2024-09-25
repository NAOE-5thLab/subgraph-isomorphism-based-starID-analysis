import os
import itertools

import numpy as np
import pandas as pd

from .yale import YaleStarCatalog
from .starDB import StarDatabase
from .utils import inter_star_angle_RADE
from utils import equatorial2vec


class PairDatabase:
    def __init__(self, D_DB, log_dir="./"):
        self.D_DB = D_DB
        self.N_DB = len(self.D_DB)
        self.D_DB_index = np.arange(self.N_DB)
        self.log_dir = log_dir
        #
        path = f"{self.log_dir}/temp.csv"
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

    def load_catalog(self):
        # pair angle
        angle_df = pd.read_csv(f"{self.log_dir}/angle_catalog.csv")
        angles_filtered = angle_df["angle"].to_numpy()
        pair_index_filtered = angle_df[["star index 1", "star index 2"]].to_numpy()
        # for index search
        df = pd.read_csv(f"{self.log_dir}/P_DB_search_index.csv")
        self.y_vector = df["y"].to_numpy()
        self.i_vector = df["i"].to_numpy()
        self.s_vector = df["s"].to_numpy()
        #
        self.angles_filtered = angles_filtered
        self.pair_index_filtered = pair_index_filtered

    def create_catalog(self, theta_max):
        # comp pair angle
        angles_filtered, pair_index_filtered = self._comp_inter_angles(theta_max)
        df = pd.DataFrame(angles_filtered, columns=["angle"])
        df["star index 1"] = pair_index_filtered[:, 0]
        df["star index 2"] = pair_index_filtered[:, 1]
        df.to_csv(f"{self.log_dir}/P_DB.csv")
        # for index search
        self.y_vector = np.array(angles_filtered)
        self.i_vector = np.argsort(self.y_vector)
        self.s_vector = self.y_vector[self.i_vector]
        df = pd.DataFrame(self.y_vector, columns=["y"])
        df["i"] = self.i_vector
        df["s"] = self.s_vector
        df.to_csv(f"{self.log_dir}/P_DB_search_index.csv")
        #
        self.angles_filtered = angles_filtered
        self.pair_index_filtered = pair_index_filtered

    def _comp_inter_angles(self, theta_max):
        # Arbitrary pairs
        index = self.D_DB_index
        p_index = np.array(list(itertools.combinations(index, 2)))
        RA1 = self.D_DB["RA [rad]"].to_numpy()[p_index[:, 0]]
        DE1 = self.D_DB["DE [rad]"].to_numpy()[p_index[:, 0]]
        RA2 = self.D_DB["RA [rad]"].to_numpy()[p_index[:, 1]]
        DE2 = self.D_DB["DE [rad]"].to_numpy()[p_index[:, 1]]
        inter_angles = inter_star_angle_RADE(RA1, DE1, RA2, DE2)
        # Pairs in theta_max
        flag_theta_max = inter_angles <= theta_max
        angles_filtered = inter_angles[flag_theta_max]
        pair_index_filtered = p_index[flag_theta_max]
        return angles_filtered, pair_index_filtered

    def get_pair_index(self):
        return self.pair_index_filtered

    def get_inter_angles(self):
        return self.angles_filtered

    def get_interval_index(self, y_lower_bound, y_upper_bound):
        lower_erase = y_lower_bound <= self.s_vector
        upper_erase = self.s_vector <= y_upper_bound
        clip = lower_erase * upper_erase
        interval_index = self.i_vector[clip]
        return interval_index

    def get_interval_pair_index(self, y_lower_bound, y_upper_bound):
        interval_index = self.get_interval_index(y_lower_bound, y_upper_bound)
        pair_index = self.pair_index_filtered[interval_index]
        return pair_index


def test():
    M_lim = 5.5
    theta_min = 1.0e-5
    theta_max = np.deg2rad(30)
    #
    # Catalog
    catalog = YaleStarCatalog()
    df_D_C = catalog.get_df()
    s_vec = equatorial2vec(catalog.get_RA(), catalog.get_DE())
    df_D_C.loc[:, ("s_X", "s_Y", "s_Z")] = s_vec
    # StarDB
    D_DB = StarDatabase(df_D_C)
    D_DB.filtering_by_visual_magnitude(M_lim)
    D_DB.filtering_by_multiple_stars(theta_min)
    D_DB.get_info()
    df_D_DB = D_DB.get_df()
    # PairDB
    P_DB = PairDatabase(df_D_DB)
    P_DB.create_catalog(theta_max)


if __name__ == "__main__":
    test()
    print("task complete")
