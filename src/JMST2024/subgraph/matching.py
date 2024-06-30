from functools import reduce
from itertools import product

import numpy as np

from .starDB import StarDatabase
from .pairDB import PairDatabase
from .utils import inter_star_angle_vec, specular_sign


class SubgraphIsomorphismBasedMatching:
    def __init__(
        self,
        D_DB: StarDatabase,
        P_DB: PairDatabase,
        eps: float,
    ) -> None:
        self.D_DB = D_DB
        self.P_DB = P_DB
        self.eps = eps

    def match_2_stars(self, s_hat_n, s_hat_m):
        theta_hat_nm = inter_star_angle_vec(s_hat_n, s_hat_m)
        candi_pairIDs = self.P_DB.get_interval_pair_index(
            theta_hat_nm - self.eps,
            theta_hat_nm + self.eps,
        )
        return candi_pairIDs

    def match_3_stars(self, s_hat_n, s_hat_m, s_hat_l):
        theta_hat_nm = inter_star_angle_vec(s_hat_n, s_hat_m)
        theta_hat_ml = inter_star_angle_vec(s_hat_m, s_hat_l)
        theta_hat_nl = inter_star_angle_vec(s_hat_n, s_hat_l)
        # filtering by each inter angle
        candi_pairIDs_nm = self.P_DB.get_interval_pair_index(
            theta_hat_nm - self.eps,
            theta_hat_nm + self.eps,
        )
        candi_pairIDs_ml = self.P_DB.get_interval_pair_index(
            theta_hat_ml - self.eps,
            theta_hat_ml + self.eps,
        )
        candi_pairIDs_nl = self.P_DB.get_interval_pair_index(
            theta_hat_nl - self.eps,
            theta_hat_nl + self.eps,
        )
        # select n
        candi_starIDs_n = np.intersect1d(
            candi_pairIDs_nm.flatten(),
            candi_pairIDs_nl.flatten(),
        )
        candi_triangleIDs_temp = []
        for candi_starID_n in candi_starIDs_n:
            starIDs_m_0 = candi_pairIDs_nm[candi_pairIDs_nm[:, 0] == candi_starID_n, 1]
            starIDs_m_1 = candi_pairIDs_nm[candi_pairIDs_nm[:, 1] == candi_starID_n, 0]
            starIDs_m_list = list(starIDs_m_0) + list(starIDs_m_1)
            starIDs_l_0 = candi_pairIDs_nl[candi_pairIDs_nl[:, 0] == candi_starID_n, 1]
            starIDs_l_1 = candi_pairIDs_nl[candi_pairIDs_nl[:, 1] == candi_starID_n, 0]
            starIDs_l_list = list(starIDs_l_0) + list(starIDs_l_1)
            star_pairIDs_ml_list = list(product(starIDs_m_list, starIDs_l_list))
            #
            for star_pairID_ml in star_pairIDs_ml_list:
                star_pairIDs_ml_sorted = np.sort(star_pairID_ml)
                matching_flag = np.all(
                    candi_pairIDs_ml == star_pairIDs_ml_sorted, axis=1
                )
                if np.any(matching_flag):
                    candi_triangleIDs_temp.append(
                        [candi_starID_n, star_pairID_ml[0], star_pairID_ml[1]]
                    )
        ### check specular sign ###
        s_vec = self.D_DB.get_s_vec()
        obs_specular_sign = specular_sign(s_hat_n, s_hat_m, s_hat_l)
        candi_triangleIDs = []
        for candi_triangleID in candi_triangleIDs_temp:
            s_vec_n = s_vec[candi_triangleID[0], :]
            s_vec_m = s_vec[candi_triangleID[1], :]
            s_vec_l = s_vec[candi_triangleID[2], :]
            if obs_specular_sign != specular_sign(s_vec_n, s_vec_m, s_vec_l):
                continue
            else:
                candi_triangleIDs.append(candi_triangleID)
        return candi_triangleIDs

    def match_4_stars(self, s_hat_n, s_hat_m, s_hat_l, s_hat_k):
        candi_triangleIDs = self.match_3_stars(s_hat_n, s_hat_m, s_hat_l)
        # matching pair stars including last star
        candi_pairIDs_ik_list = []
        for s_hat in [s_hat_n, s_hat_m, s_hat_l]:
            theta_hat_ik = inter_star_angle_vec(s_hat, s_hat_k)
            candi_pairIDs_ik = self.P_DB.get_interval_pair_index(
                theta_hat_ik - self.eps,
                theta_hat_ik + self.eps,
            )
            candi_pairIDs_ik_list.append(candi_pairIDs_ik)
        # determine candidate of last star in each set of matched stars except last
        candi_pyramidIDs = []
        for candi_triangleID in candi_triangleIDs:
            candi_pyramidIDs_k_from_i_list = []
            for i, candi_pairIDs_ik in enumerate(candi_pairIDs_ik_list):
                candi_starID_i = candi_triangleID[i]
                # determine candidate of last star from inter angle between last and i-th star
                match_i0 = candi_pairIDs_ik[:, 0] == candi_starID_i
                candi_pyramidIDs_k_from_i0 = candi_pairIDs_ik[match_i0, 1]
                match_i1 = candi_pairIDs_ik[:, 1] == candi_starID_i
                candi_pyramidIDs_k_from_i1 = candi_pairIDs_ik[match_i1, 0]
                candi_pyramidIDs_k_from_i = np.union1d(
                    candi_pyramidIDs_k_from_i0, candi_pyramidIDs_k_from_i1
                )
                candi_pyramidIDs_k_from_i_list.append(candi_pyramidIDs_k_from_i)
            candi_pyramidIDs_k = reduce(
                np.intersect1d, tuple(candi_pyramidIDs_k_from_i_list)
            )
            for candi_pyramidID_k in candi_pyramidIDs_k:
                candi_pyramidIDs.append(list(candi_triangleID) + [candi_pyramidID_k])
        return candi_pyramidIDs

    def match_p_stars(self, s_hat_list):
        p = len(s_hat_list)
        if p < 2:
            return [], range(p)
        elif p == 2:
            candi_pairIDs = self.match_2_stars(*s_hat_list)
            return candi_pairIDs, range(p)
        elif p == 3:
            candi_triangleIDs = self.match_3_stars(*s_hat_list)
            return candi_triangleIDs, range(p)
        elif p == 4:
            candi_pyramidIDs = self.match_4_stars(*s_hat_list)
            return candi_pyramidIDs, range(p)
        else:
            for i, j, k, l in self.pattern_shifting_4(p):
                s_hat_list_ = [s_hat_list[n] for n in [i, j, k, l]]
                candi_pyramidIDs = self.match_4_stars(*s_hat_list_)
                if len(candi_pyramidIDs) == 1:
                    return candi_pyramidIDs, [i, j, k, l]
            return candi_pyramidIDs, [i, j, k, l]

    def pattern_shifting_3(self, N):
        for di in range(1, N - 2 + 1):
            for dj in range(1, N - di - 1 + 1):
                for i in range(N - di - dj):
                    j = i + di
                    k = j + dj
                    yield i, j, k

    def pattern_shifting_4(self, N):
        for di in range(1, N - 3 + 1):
            for dj in range(1, N - di - 2 + 1):
                for dk in range(1, N - di - dj - 1 + 1):
                    for i in range(N - di - dj - dk):
                        j = i + di
                        k = j + dj
                        l = k + dk
                        yield i, j, k, l
