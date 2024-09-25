import itertools
import numpy as np

from database.yale import YaleStarCatalog
from utils import inter_star_angle_RADE


class StarCatalog:
    def __init__(self, ID, RA, DE):
        self.ID = ID
        self.RA = RA
        self.DE = DE
        self.N = len(self.ID)
        self.index = np.arange(self.N)
        #
        self.filter = np.array([True]*len(ID))

    def get_ID(self):
        return self.ID[self.filter]

    def get_RA(self):
        return self.RA[self.filter]

    def get_DE(self):
        return self.DE[self.filter]

    def filtering_by_visual_magnitude(self, Mv, Mv_max):
        Mv_filter = Mv <= Mv_max
        self.filter = self.filter * Mv_filter

    def filtering_by_variable_stars(self, Var):
        Var_filter = Var == ''
        self.filter = self.filter * Var_filter

    def filtering_by_multiple_stars(self, multi_angle):
        index = self.index[self.filter]
        pair_index = np.array(list(itertools.combinations(index, 2)))
        inter_angles = inter_star_angle_RADE(
            self.RA[pair_index[:, 0]], self.DE[pair_index[:, 0]],
            self.RA[pair_index[:, 1]], self.DE[pair_index[:, 1]])
        multi_index = np.unique(
            pair_index[inter_angles < multi_angle].flatten())
        multi_filter = np.array(
            [False if i in multi_index else True for i in self.index])
        self.filter = self.filter * multi_filter


def test():
    Mv_max = 5.0
    multi_angle = 1.0e-5
    #
    yale = YaleStarCatalog()
    star_catalog = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    star_catalog.filtering_by_visual_magnitude(yale.get_Vmag(), Mv_max)
    star_catalog.filtering_by_variable_stars(yale.get_VarID())
    star_catalog.filtering_by_multiple_stars(multi_angle)
    ID = star_catalog.get_ID()
    print(ID)


if __name__ == '__main__':
    test()
    print('task completed')
