from functools import reduce
from itertools import product, combinations

import numpy as np

from .starDB import StarDatabase
from .pairDB import PairDatabase
from .utils import inter_star_angle_vec, specular_sign


class SubgraphIsomorphismBasedMatching:
    """
    Class to perform subgraph isomorphism-based matching of stars.
    """

    def __init__(
        self,
        D_DB: StarDatabase,
        P_DB: PairDatabase,
        eps: float,
        np_random: np.random.Generator = None,
    ) -> None:
        """
        Initialize the SubgraphIsomorphismBasedMatching class.

        Args:
            D_DB (StarDatabase): Star database.
            P_DB (PairDatabase): Pair database.
            eps (float): Tolerance for matching angles.
        """
        self.D_DB = D_DB
        self.P_DB = P_DB
        self.eps = eps
        self.np_random = np.random.default_rng() if np_random is None else np_random

    def match_stars(self, s_hat_list: list, with_check: bool = False):
        """
        Match stars based on the given list of star vectors.

        Args:
            s_hat_list (list): List of star vectors.
            with_check (bool): Whether to check the matching result.

        Returns:
            tuple: Candidate star combinations, observation IDs, and matching info.
        """
        m2to3 = self.match_3_stars_with_matched_2_stars
        mp1top = self.match_p_stars_with_matched_pminus1_stars
        #
        info = ""
        if len(s_hat_list) == 0:
            return [], [], "Given no star"
        for obs_IDs, pattern in self.pattern_shuffling(s_hat_list):
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
            info += f"Matching number of {obs_IDs} : {len(candi_pcombIDs)}\n"
            # update
            if (len(candi_pcombIDs) <= 1) and (p >= 2):
                break
            candi_pminus1combIDs = candi_pcombIDs
        return candi_pcombIDs, obs_IDs, info

    def match_2_stars(self, s_hat_list):
        """
        Match two stars based on their angular separation.

        Args:
            s_hat_list (list): List of two star vectors.

        Returns:
            np.ndarray: Candidate pair IDs.
        """
        assert len(s_hat_list) == 2
        theta_hat_nm = inter_star_angle_vec(*s_hat_list)
        candi_pairIDs = self.P_DB.get_interval_pair_index(
            theta_hat_nm - self.eps,
            theta_hat_nm + self.eps,
        )
        return candi_pairIDs

    def match_3_stars_with_matched_2_stars(self, s_hat_list, candi_pairIDs_nm):
        """
        Match three stars based on their angular separations and matched pairs.

        Args:
            s_hat_list (list): List of three star vectors.
            candi_pairIDs_nm (np.ndarray): Candidate pair IDs for the first two stars.

        Returns:
            list: Candidate triangle IDs.
        """
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
        # Check specular sign
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
        """
        Match three stars based on their angular separations.

        Args:
            s_hat_list (list): List of three star vectors.

        Returns:
            list: Candidate triangle IDs.
        """
        assert len(s_hat_list) == 3
        candi_pairIDs_nm = self.match_2_stars(s_hat_list[:2])
        candi_triangleIDs = self.match_3_stars_with_matched_2_stars(
            s_hat_list, candi_pairIDs_nm
        )
        return candi_triangleIDs

    def match_p_stars_with_matched_pminus1_stars(
        self, s_hat_list, candi_pminus1combIDs
    ):
        """
        Match p stars based on their angular separations and matched (p-1) stars.

        Args:
            s_hat_list (list): List of p star vectors.
            candi_pminus1combIDs (list): Candidate combinations of (p-1) stars.

        Returns:
            list: Candidate combinations of p stars.
        """
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
        """
        Match p stars based on their angular separations.

        Args:
            s_hat_list (list): List of p star vectors.

        Returns:
            list: Candidate combinations of p stars.
        """
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

    def pattern_shuffling(self, pattern: list):
        """
        Shuffle the pattern list.

        Args:
            pattern (list): List of patterns.

        Yields:
            tuple: Selected index and selected pattern.
        """
        N = len(pattern)
        selected_index = []
        selected_pattern = []
        remained_index = list(range(N))
        for n in range(N):
            i = self.np_random.choice(remained_index)
            selected_index.append(i)
            selected_pattern.append(pattern[i])
            remained_index.remove(i)
            yield selected_index, selected_pattern

    def check_matching_result(self, candi_pcombIDs, s_hat_list):
        """
        Check the matching result.

        Args:
            candi_pcombIDs (list): Candidate combinations of p stars.
            s_hat_list (list): List of star vectors.

        Returns:
            bool: True if the matching result is correct, False otherwise.
        """
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
