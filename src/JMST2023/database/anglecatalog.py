import os
import itertools

import numpy as np
import pandas as pd

from database.yale import YaleStarCatalog
from database.catalog import StarCatalog
from database.indexsearch import IndexSearch
from utils.inter_star_angle import inter_star_angle_RADE
from utils.unit_convert import deg2rad

dir = os.path.dirname(__file__)


class InterAngleCatalog:
    def __init__(self, ID, RA, DE, log_dir='./'):
        self.ID = ID
        self.RA = RA
        self.DE = DE
        self.N = len(self.ID)
        self.index = range(self.N)
        #
        self.log_dir = log_dir

    def create_catalog(self, FOV=deg2rad(30), use_log=False):
        if use_log:
            angle_df = pd.read_csv(f'{self.log_dir}/angle_catalog.csv')
            angles_filtered = angle_df['angle'].to_numpy()
            pair_index_filtered = angle_df[[
                'star index 1', 'star index 2']].to_numpy()
            #
            index_search = IndexSearch()
            index_search.load_vector(self.log_dir, 'angle_catalog_search_index.csv')
        else:
            # Arbitrary pairs
            pair_index = np.array(
                list(itertools.combinations(self.index, 2)))
            inter_angles = inter_star_angle_RADE(
                self.RA[pair_index[:, 0]], self.DE[pair_index[:, 0]],
                self.RA[pair_index[:, 1]], self.DE[pair_index[:, 1]])
            # Pairs in FOV
            filter_FOV = inter_angles < FOV
            angles_filtered = inter_angles[filter_FOV]
            pair_index_filtered = pair_index[filter_FOV]
            angle_df = pd.DataFrame(angles_filtered, columns=['angle'])
            angle_df['star index 1'] = pair_index_filtered[:, 0]
            angle_df['star index 2'] = pair_index_filtered[:, 1]
            angle_df.to_csv(f'{self.log_dir}/angle_catalog.csv')
            #
            index_search = IndexSearch(y_vector=angles_filtered)
            index_search.save_vector(self.log_dir, 'angle_catalog_search_index.csv')
        #
        self.angles_filtered = angles_filtered
        self.pair_index_filtered = pair_index_filtered
        self.index_search = index_search

    def get_pair_index(self):
        return self.pair_index_filtered

    def get_inter_angles(self):
        return self.angles_filtered

    def get_interval_index(self, y_lower_bound, y_upper_bound):
        interval_index = self.index_search.get_interval_index(
            y_lower_bound, y_upper_bound)
        return interval_index

    def get_interval_pair_index(self, y_lower_bound, y_upper_bound):
        interval_index = self.index_search.get_interval_index(
            y_lower_bound, y_upper_bound)
        pair_index = self.pair_index_filtered[interval_index]
        return pair_index


def test():
    Mv_max = 5.5
    multi_angle = 1.0e-5
    FOV = deg2rad(30)
    #
    yale = YaleStarCatalog()
    star_ctlg = StarCatalog(yale.get_HR(), yale.get_RA(), yale.get_DE())
    star_ctlg.filtering_by_visual_magnitude(yale.get_Vmag(), Mv_max)
    star_ctlg.filtering_by_variable_stars(yale.get_VarID())
    star_ctlg.filtering_by_multiple_stars(multi_angle)
    angle_ctlg = InterAngleCatalog(
        star_ctlg.get_ID(), star_ctlg.get_RA(), star_ctlg.get_DE())
    angle_ctlg.create_catalog(FOV=FOV, use_log=False)


if __name__ == '__main__':
    test()
    print('task complete')
