import itertools
import random
from functools import reduce

import numpy as np

from utils import specular_sign, equatorial2vec, inter_star_angle_vec


def search_pairIDs_of_matched_theta(theta, epsillon, pairstar_db):
    # matching pair graph
    candi_pairIDs = pairstar_db.get_pairIDs_in_certain_interval(
        theta - epsillon,
        theta + epsillon)
    return candi_pairIDs


def matching_pair(s_hat_list, epsillon, pairstar_db):
    # calc inter angle of observed stars
    s_hat_i = s_hat_list[0]
    s_hat_j = s_hat_list[1]
    theta_hat_ij = inter_star_angle_vec(s_hat_i, s_hat_j)
    # matching pair graph
    candi_pairIDs = search_pairIDs_of_matched_theta(
        theta_hat_ij, epsillon, pairstar_db)
    return candi_pairIDs


def matching_triangle(s_hat_list, epsillon, pairstar_db):
    RA = pairstar_db.get_RA()
    DE = pairstar_db.get_DE()
    s_hat_i = s_hat_list[0]
    s_hat_j = s_hat_list[1]
    s_hat_k = s_hat_list[2]
    theta_hat_ij = inter_star_angle_vec(s_hat_i, s_hat_j)
    theta_hat_ik = inter_star_angle_vec(s_hat_i, s_hat_k)
    theta_hat_jk = inter_star_angle_vec(s_hat_j, s_hat_k)
    ### select triangle ###
    # filtering by each inter angle
    candi_pairIDs_ij = search_pairIDs_of_matched_theta(
        theta_hat_ij, epsillon, pairstar_db)
    candi_pairIDs_ik = search_pairIDs_of_matched_theta(
        theta_hat_ik, epsillon, pairstar_db)
    candi_pairIDs_jk = search_pairIDs_of_matched_theta(
        theta_hat_jk, epsillon, pairstar_db)
    # select i
    candi_starIDs_i = np.intersect1d(
        candi_pairIDs_ij.flatten(), candi_pairIDs_ik.flatten())
    candi_triangleIDs_temp = []
    for candi_starID_i in candi_starIDs_i:
        starIDs_j_0 = candi_pairIDs_ij[candi_pairIDs_ij[:, 0]
                                       == candi_starID_i, 1]
        starIDs_j_1 = candi_pairIDs_ij[candi_pairIDs_ij[:, 1]
                                       == candi_starID_i, 0]
        starIDs_k_0 = candi_pairIDs_ik[candi_pairIDs_ik[:, 0]
                                       == candi_starID_i, 1]
        starIDs_k_1 = candi_pairIDs_ik[candi_pairIDs_ik[:, 1]
                                       == candi_starID_i, 0]
        starIDs_k_list = list(starIDs_k_0) + list(starIDs_k_1)
        starIDs_j_list = list(starIDs_j_0) + list(starIDs_j_1)
        star_pairIDs_jk_list = list(itertools.product(
            starIDs_j_list, starIDs_k_list))
        #
        for star_pairID_jk in star_pairIDs_jk_list:
            star_pairIDs_jk_sorted = np.sort(star_pairID_jk)
            matching_flag = np.all(candi_pairIDs_jk ==
                                   star_pairIDs_jk_sorted, axis=1)
            if np.any(matching_flag):
                candi_triangleIDs_temp.append(
                    [candi_starID_i, star_pairID_jk[0], star_pairID_jk[1]])    #
    ### check specular sign ###
    obs_specular_sign = specular_sign(s_hat_i, s_hat_j, s_hat_k)
    candi_triangleIDs = []
    for candi_triangleID in candi_triangleIDs_temp:
        RA_i = RA[candi_triangleID[0]]
        DE_i = DE[candi_triangleID[0]]
        RADE_vec_i = equatorial2vec(RA_i, DE_i)
        RA_j = RA[candi_triangleID[1]]
        DE_j = DE[candi_triangleID[1]]
        RADE_vec_j = equatorial2vec(RA_j, DE_j)
        RA_k = RA[candi_triangleID[2]]
        DE_k = DE[candi_triangleID[2]]
        RADE_vec_k = equatorial2vec(RA_k, DE_k)
        #
        catalog_specular_sign = specular_sign(
            RADE_vec_i, RADE_vec_j, RADE_vec_k)
        if obs_specular_sign != catalog_specular_sign:
            continue
        else:
            candi_triangleIDs.append(candi_triangleID)
    # candi_triangleIDs = np.array(candi_triangleIDs)
    return candi_triangleIDs


def matching_subgraph(n_obs, s_hat_list, epsillon, pairstar_db):
    if n_obs > 3:
        s_hat_list_except_last = s_hat_list[:-1]
        s_hat_last = s_hat_list[-1]
        # matching stars except last
        candi_indexes_except_last = matching_subgraph(
            n_obs - 1, s_hat_list_except_last, epsillon, pairstar_db)
        # matching pair stars including last star
        candi_star_indexes_ilast_list = []
        for i in range(n_obs-1):
            theta_hat_ilast = inter_star_angle_vec(
                s_hat_list_except_last[i], s_hat_last)
            candi_star_indexes_ilast = search_pairIDs_of_matched_theta(
                theta_hat_ilast, epsillon, pairstar_db)
            candi_star_indexes_ilast_list.append(candi_star_indexes_ilast)
        # determine candidate of last star in each set of matched stars except last
        candi_indexes = []
        for candi_index_except_last in candi_indexes_except_last:
            candi_star_indexes_last_from_i_list = []
            for i in range(n_obs-1):
                # define by other name
                candi_star_indexes_ilast = candi_star_indexes_ilast_list[i]
                candi_star_i = candi_index_except_last[i]
                # determine candidate of last star from inter angle between last and i-th star
                candi_star_indexes_last_from_i0 = candi_star_indexes_ilast[
                    candi_star_indexes_ilast[:, 0] == candi_star_i, 1]
                candi_star_indexes_last_from_i1 = candi_star_indexes_ilast[
                    candi_star_indexes_ilast[:, 1] == candi_star_i, 0]
                candi_star_indexes_last_from_i = np.union1d(
                    candi_star_indexes_last_from_i0, candi_star_indexes_last_from_i1)
                candi_star_indexes_last_from_i_list.append(
                    candi_star_indexes_last_from_i)

            candi_star_indexes_last = reduce(
                np.intersect1d, tuple(candi_star_indexes_last_from_i_list))
            #
            for candi_star_index_last in candi_star_indexes_last:
                candi_indexes.append(
                    list(candi_index_except_last) + [candi_star_index_last])
        return candi_indexes
    elif n_obs > 2:
        return matching_triangle(s_hat_list, epsillon, pairstar_db)
    else:
        return matching_pair(s_hat_list, epsillon, pairstar_db)


# def matching_pyramid(
#         theta_hat_ij, theta_hat_ik, theta_hat_jk,
#         theta_hat_il, theta_hat_jl, theta_hat_kl,
#         s_hat_i, s_hat_j, s_hat_k, s_hat_l,
#         epsillon, pairstar_db):
#     ### select pyramid ###
#     candi_triangle_indexes = matching_triangle(
#         theta_hat_ij, theta_hat_ik, theta_hat_jk,
#         s_hat_i, s_hat_j, s_hat_k,
#         epsillon, pairstar_db)
#     # filtering by each inter angle
#     star_indexes_il = matching_pair(theta_hat_il, epsillon, pairstar_db)
#     star_indexes_jl = matching_pair(theta_hat_jl, epsillon, pairstar_db)
#     star_indexes_kl = matching_pair(theta_hat_kl, epsillon, pairstar_db)
#     #
#     candi_pyramid_indexes = []
#     for candi_triangle_index in candi_triangle_indexes:
#         #
#         candi_star_i = candi_triangle_index[0]
#         candi_star_j = candi_triangle_index[1]
#         candi_star_k = candi_triangle_index[2]
#         #
#         candi_star_indexes_l_from_i0 = star_indexes_il[star_indexes_il[:, 0]
#                                                        == candi_star_i, 1]
#         candi_star_indexes_l_from_i1 = star_indexes_il[star_indexes_il[:, 1]
#                                                        == candi_star_i, 0]
#         candi_star_indexes_l_from_i = np.union1d(
#             candi_star_indexes_l_from_i0, candi_star_indexes_l_from_i1)
#         candi_star_indexes_l_from_j0 = star_indexes_jl[star_indexes_jl[:, 0]
#                                                        == candi_star_j, 1]
#         candi_star_indexes_l_from_j1 = star_indexes_jl[star_indexes_jl[:, 1]
#                                                        == candi_star_j, 0]
#         candi_star_indexes_l_from_j = np.union1d(
#             candi_star_indexes_l_from_j0, candi_star_indexes_l_from_j1)
#         candi_star_indexes_l_from_k0 = star_indexes_kl[star_indexes_kl[:, 0]
#                                                        == candi_star_k, 1]
#         candi_star_indexes_l_from_k1 = star_indexes_kl[star_indexes_kl[:, 1]
#                                                        == candi_star_k, 0]
#         candi_star_indexes_l_from_k = np.union1d(
#             candi_star_indexes_l_from_k0, candi_star_indexes_l_from_k1)
#         #
#         candi_star_indexes_l = reduce(
#             np.intersect1d,
#             (candi_star_indexes_l_from_i,
#              candi_star_indexes_l_from_j, candi_star_indexes_l_from_k)
#         )
#         #
#         for candi_star_index_l in candi_star_indexes_l:
#             candi_pyramid_indexes.append([
#                 candi_triangle_index[0],
#                 candi_triangle_index[1],
#                 candi_triangle_index[2],
#                 candi_star_index_l])
#     return np.array(candi_pyramid_indexes)


# def calc_Theta_pair_hat(n_obs, s_hat_list):
#     n_obs = len(s_hat_list)
#     pair_indexes = list(itertools.combinations(range(n_obs), 2))
#     theta_hat_list = []
#     for pair_index in pair_indexes:
#         s1_hat = s_hat_list[pair_index[0]]
#         s2_hat = s_hat_list[pair_index[1]]
#         theta_hat = inter_star_angle_vec(s1_hat, s2_hat)
#         theta_hat_list.append(theta_hat)
#     return theta_hat_list
