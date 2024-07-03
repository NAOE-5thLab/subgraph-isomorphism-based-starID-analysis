from functools import reduce
from itertools import product, combinations

import numpy as np

from .starDB import StarDatabase
from .pairDB import PairDatabase
from .utils import inter_star_angle_vec, specular_sign


class SubgraphIsomorphismBasedMatching:
    def __init__(self, D_DB: StarDatabase, P_DB: PairDatabase, eps: float) -> None:
        self.D_DB = D_DB
        self.P_DB = P_DB
        self.eps = eps

    def match_stars(
        self,
        s_hat_list: list,
        np_random: np.random.Generator,
        with_check: bool = False,
    ):
        m2to3 = self.match_3_stars_with_matched_2_stars
        mp1top = self.match_p_stars_with_matched_pminus1_stars
        #
        info = None
        if len(s_hat_list) == 0:
            return [], [], {"N_candi": []}
        for obs_IDs, pattern in self.pattern_shuffling(s_hat_list, np_random):
            p = len(pattern)
            if p == 1:
                candi_pcombIDs = [[i] for i in range(len(self.D_DB.D_DB))]
            elif p == 2:
                candi_pcombIDs = self.match_2_stars(pattern)
            elif p == 3:
                candi_pcombIDs = m2to3(pattern, candi_pminus1combIDs)
            else:
                candi_pcombIDs = mp1top(pattern, candi_pminus1combIDs)
            # check
            if with_check:
                assert self.check_matching_result(candi_pcombIDs, pattern)
            # info
            if info is None:
                info = {}
                info["N_candi"] = []
            info["N_candi"].append(len(candi_pcombIDs))
            # update
            if (len(candi_pcombIDs) <= 1) and (p >= 2):
                break
            candi_pminus1combIDs = candi_pcombIDs
        return candi_pcombIDs, obs_IDs, info

    def match_2_stars(self, s_hat_list):
        assert len(s_hat_list) == 2
        theta_hat_nm = inter_star_angle_vec(*s_hat_list)
        candi_pairIDs = self.P_DB.get_interval_pair_index(
            theta_hat_nm - self.eps,
            theta_hat_nm + self.eps,
        )
        return candi_pairIDs

    def match_3_stars_with_matched_2_stars(self, s_hat_list, candi_pairIDs_nm):
        assert len(s_hat_list) == 3
        s_hat_n, s_hat_m, s_hat_l = s_hat_list
        # compute inter angle
        theta_hat_ml = inter_star_angle_vec(s_hat_m, s_hat_l)
        theta_hat_nl = inter_star_angle_vec(s_hat_n, s_hat_l)
        # filtering by each inter angle
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

    def match_3_stars(self, s_hat_list):
        assert len(s_hat_list) == 3
        candi_pairIDs_nm = self.match_2_stars(s_hat_list[:2])
        candi_triangleIDs = self.match_3_stars_with_matched_2_stars(
            s_hat_list, candi_pairIDs_nm
        )
        return candi_triangleIDs

    def match_p_stars_with_matched_pminus1_stars(
        self, s_hat_list, candi_pminus1combIDs
    ):
        assert len(s_hat_list) - 1 == len(candi_pminus1combIDs[0])
        s_hat_list_except_p = s_hat_list[:-1]
        s_hat_p = s_hat_list[-1]
        # matching pair stars including last star
        candi_pairIDs_ip_list = []
        for s_hat in s_hat_list_except_p:
            theta_hat_ip = inter_star_angle_vec(s_hat, s_hat_p)
            candi_pairIDs_ip = self.P_DB.get_interval_pair_index(
                theta_hat_ip - self.eps,
                theta_hat_ip + self.eps,
            )
            candi_pairIDs_ip_list.append(candi_pairIDs_ip)
        # determine candidate of last star in each set of matched stars except last
        candi_pcombIDs = []
        for candi_pminus1combID in candi_pminus1combIDs:
            candi_pcombIDs_p_from_i_list = []
            for i, candi_pairIDs_ip in enumerate(candi_pairIDs_ip_list):
                candi_starID_i = candi_pminus1combID[i]
                # determine candidate of last star from inter angle between last and i-th star
                match_i0 = candi_pairIDs_ip[:, 0] == candi_starID_i
                candi_pcombIDs_p_from_i0 = candi_pairIDs_ip[match_i0, 1]
                match_i1 = candi_pairIDs_ip[:, 1] == candi_starID_i
                candi_pcombIDs_p_from_i1 = candi_pairIDs_ip[match_i1, 0]
                candi_pcombIDs_p_from_i = np.union1d(
                    candi_pcombIDs_p_from_i0, candi_pcombIDs_p_from_i1
                )
                candi_pcombIDs_p_from_i_list.append(candi_pcombIDs_p_from_i)
            candi_pcombIDs_p = reduce(
                np.intersect1d, tuple(candi_pcombIDs_p_from_i_list)
            )
            for candi_pcombID_p in candi_pcombIDs_p:
                candi_pcombIDs.append(list(candi_pminus1combID) + [candi_pcombID_p])
        return candi_pcombIDs

    def match_p_stars(self, s_hat_list):
        assert len(s_hat_list) > 3
        match_2to3 = self.match_3_stars_with_matched_2_stars
        match_p1top = self.match_p_stars_with_matched_pminus1_stars
        #
        candi_pminus1combIDs = self.match_2_stars(s_hat_list[:2])
        candi_pminus1combIDs = match_2to3(s_hat_list[:3], candi_pminus1combIDs)
        for pi in range(4, len(s_hat_list)):
            candi_pcombIDs = match_p1top(s_hat_list[:pi], candi_pminus1combIDs)
            candi_pminus1combIDs = candi_pcombIDs
        return candi_pcombIDs

    def pattern_shuffling(self, pattern: list, np_random: np.random.Generator):
        N = len(pattern)
        selected_index = []
        selected_pattern = []
        remained_index = list(range(N))
        for n in range(N):
            i = np_random.choice(remained_index)
            selected_index.append(i)
            selected_pattern.append(pattern[i])
            remained_index.remove(i)
            yield selected_index, selected_pattern

    def check_matching_result(self, candi_pcombIDs, s_hat_list):
        if len(s_hat_list) <= 1:
            return True
        #
        obs_angles = []
        for s1_hat, s2_hat in list(combinations(s_hat_list, 2)):
            obs_angles.append(inter_star_angle_vec(s1_hat, s2_hat))
        for candi_pcombID in candi_pcombIDs:
            s_candi_list = self.D_DB.get_s_vec()[candi_pcombID, :]
            for i, (s1, s2) in enumerate(list(combinations(s_candi_list, 2))):
                condi_angle = inter_star_angle_vec(s1, s2)
                ok = np.abs(obs_angles[i] - condi_angle) < self.eps
                if not ok:
                    return False
        return True

    # def pattern_shifting_3(self, N):
    #     for di in range(1, N - 2 + 1):
    #         for dj in range(1, N - di - 1 + 1):
    #             for i in range(N - di - dj):
    #                 j = i + di
    #                 k = j + dj
    #                 yield i, j, k

    # def pattern_shifting_4(self, N):
    #     for di in range(1, N - 3 + 1):
    #         for dj in range(1, N - di - 2 + 1):
    #             for dk in range(1, N - di - dj - 1 + 1):
    #                 for i in range(N - di - dj - dk):
    #                     j = i + di
    #                     k = j + dj
    #                     l = k + dk
    #                     yield i, j, k, l

    # def match_4_stars(self, s_hat_n, s_hat_m, s_hat_l, s_hat_k):
    #     candi_triangleIDs = self.match_3_stars(s_hat_n, s_hat_m, s_hat_l)
    #     # matching pair stars including last star
    #     candi_pairIDs_ik_list = []
    #     for s_hat in [s_hat_n, s_hat_m, s_hat_l]:
    #         theta_hat_ik = inter_star_angle_vec(s_hat, s_hat_k)
    #         candi_pairIDs_ik = self.P_DB.get_interval_pair_index(
    #             theta_hat_ik - self.eps,
    #             theta_hat_ik + self.eps,
    #         )
    #         candi_pairIDs_ik_list.append(candi_pairIDs_ik)
    #     # determine candidate of last star in each set of matched stars except last
    #     candi_pyramidIDs = []
    #     for candi_triangleID in candi_triangleIDs:
    #         candi_pyramidIDs_k_from_i_list = []
    #         for i, candi_pairIDs_ik in enumerate(candi_pairIDs_ik_list):
    #             candi_starID_i = candi_triangleID[i]
    #             # determine candidate of last star from inter angle between last and i-th star
    #             match_i0 = candi_pairIDs_ik[:, 0] == candi_starID_i
    #             candi_pyramidIDs_k_from_i0 = candi_pairIDs_ik[match_i0, 1]
    #             match_i1 = candi_pairIDs_ik[:, 1] == candi_starID_i
    #             candi_pyramidIDs_k_from_i1 = candi_pairIDs_ik[match_i1, 0]
    #             candi_pyramidIDs_k_from_i = np.union1d(
    #                 candi_pyramidIDs_k_from_i0, candi_pyramidIDs_k_from_i1
    #             )
    #             candi_pyramidIDs_k_from_i_list.append(candi_pyramidIDs_k_from_i)
    #         candi_pyramidIDs_k = reduce(
    #             np.intersect1d, tuple(candi_pyramidIDs_k_from_i_list)
    #         )
    #         for candi_pyramidID_k in candi_pyramidIDs_k:
    #             candi_pyramidIDs.append(list(candi_triangleID) + [candi_pyramidID_k])
    #     return candi_pyramidIDs

    # def match_p_stars(self, s_hat_list, obs_HR):
    #     p = len(s_hat_list)
    #     if p < 2:
    #         return [], range(p)
    #     elif p == 2:
    #         candi_pairIDs = self.match_2_stars(*s_hat_list)
    #         return candi_pairIDs, range(p)
    #     elif p == 3:
    #         candi_triangleIDs = self.match_3_stars(*s_hat_list)
    #         return candi_triangleIDs, range(p)
    #     elif p == 4:
    #         candi_pyramidIDs = self.match_4_stars(*s_hat_list)
    #         return candi_pyramidIDs, range(p)
    #     else:
    #         for i, j, k, l in self.pattern_shifting_4(p):
    #             s_hat_list_ = [s_hat_list[n] for n in [i, j, k, l]]
    #             candi_pyramidIDs = self.match_4_stars(*s_hat_list_)
    #             print(len(candi_pyramidIDs))
    #             if len(candi_pyramidIDs) == 1:
    #                 return candi_pyramidIDs, [i, j, k, l]
    #             else:
    #                 a = np.array(
    #                     [
    #                         inter_star_angle_vec(s_hat_list_[0], s_hat_list_[1]),
    #                         inter_star_angle_vec(s_hat_list_[1], s_hat_list_[2]),
    #                         inter_star_angle_vec(s_hat_list_[2], s_hat_list_[3]),
    #                         inter_star_angle_vec(s_hat_list_[3], s_hat_list_[0]),
    #                     ]
    #                 )
    #                 print(self.eps)
    #                 for candi_ID in candi_pyramidIDs:
    #                     s = self.D_DB.get_s_vec()
    #                     b = np.array(
    #                         [
    #                             inter_star_angle_vec(s[candi_ID[0]], s[candi_ID[1]]),
    #                             inter_star_angle_vec(s[candi_ID[1]], s[candi_ID[2]]),
    #                             inter_star_angle_vec(s[candi_ID[2]], s[candi_ID[3]]),
    #                             inter_star_angle_vec(s[candi_ID[3]], s[candi_ID[0]]),
    #                         ]
    #                     )
    #                     print(np.abs(a - b))
    #                     print("")
    #         return candi_pyramidIDs, [i, j, k, l]
