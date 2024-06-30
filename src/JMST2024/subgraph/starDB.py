import itertools
import numpy as np
import pandas as pd

from .yale import YaleStarCatalog
from .utils import inter_star_angle_RADE
from utils import equatorial2vec


class StarDatabase:
    def __init__(self, D_C: pd.DataFrame):
        self.D_C = D_C
        self.N_C = len(self.D_C)
        self.D_C_index = np.arange(self.N_C)
        #
        self.DB_flag = np.array([True] * self.N_C)
        self.D_DB = self.D_C

    def filtering_by_visual_magnitude(self, M_lim, lb=None):
        if lb is not None:
            self.DB_flag *= lb < self.D_C["Vmag [mag]"].to_numpy()
        self.DB_flag *= self.D_C["Vmag [mag]"].to_numpy() < M_lim
        self.D_DB = self.D_C[self.DB_flag]

    def filtering_by_multiple_stars(self, theta_max):
        index = self.D_C_index[self.DB_flag]
        p_index = np.array(list(itertools.combinations(index, 2)))
        RA1 = self.D_C["RA [rad]"].to_numpy()[p_index[:, 0]]
        DE1 = self.D_C["DE [rad]"].to_numpy()[p_index[:, 0]]
        RA2 = self.D_C["RA [rad]"].to_numpy()[p_index[:, 1]]
        DE2 = self.D_C["DE [rad]"].to_numpy()[p_index[:, 1]]
        inter_angles = inter_star_angle_RADE(RA1, DE1, RA2, DE2)
        multi_index = np.unique(p_index[inter_angles < theta_max].flatten())
        self.DB_flag[multi_index] = False
        self.D_DB = self.D_C[self.DB_flag]

    def get_info(self):
        print("----- Star Database -----")
        Vmag = self.get_Vmag()
        print(f"the number of stars : {len(self.get_HR())}")
        print(f"the range of magnitude : [{min(Vmag)}, {max(Vmag)}]")

    def get_df(self):
        return self.D_DB[
            ["HR", "RA [rad]", "DE [rad]", "Vmag [mag]", "s_X", "s_Y", "s_Z"]
        ]

    def get_D(self):
        return self.get_df().to_numpy()

    def get_HR(self):
        return self.D_DB["HR"].to_numpy()

    def get_RA(self):
        return self.D_DB["RA [rad]"].to_numpy()

    def get_DE(self):
        return self.D_DB["DE [rad]"].to_numpy()

    def get_s_vec(self):
        return self.D_DB[["s_X", "s_Y", "s_Z"]].to_numpy()

    def get_Vmag(self):
        return self.D_DB["Vmag [mag]"].to_numpy()


def test():
    M_lim = 5.0
    theta_min = 1.0e-5
    #
    catalog = YaleStarCatalog()
    df_D_C = catalog.get_df()
    df_D_C[["s_X", "s_Y", "s_Z"]] = equatorial2vec(catalog.get_RA(), catalog.get_DE())
    # StarDB
    D_DB = StarDatabase(df_D_C)
    D_DB.filtering_by_visual_magnitude(M_lim)
    D_DB.filtering_by_multiple_stars(theta_min)
    ID = D_DB.get_HR()
    print(ID)


if __name__ == "__main__":
    test()
    print("task completed")
