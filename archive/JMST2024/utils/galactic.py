import numpy as np


alpha_NGP = np.deg2rad(192.85948)
delta_NGP = np.deg2rad(27.12825)
l_NGP = np.deg2rad(122.93192)


def equatorial2galactic(alpha, delta):
    sin_b_1 = np.sin(delta_NGP) * np.sin(delta)
    sin_b_2 = np.cos(delta_NGP) * np.cos(delta) * np.cos(alpha - alpha_NGP)
    sin_b = sin_b_1 + sin_b_2
    cos_b_sin_l_NGP_l = np.cos(delta) * np.sin(alpha - alpha_NGP)
    cos_b_cos_l_NGP_l_1 = np.cos(delta_NGP) * np.sin(delta)
    cos_b_cos_l_NGP_l_2 = np.sin(delta_NGP) * np.cos(delta) * np.cos(alpha - alpha_NGP)
    cos_b_cos_l_NGP_l = cos_b_cos_l_NGP_l_1 - cos_b_cos_l_NGP_l_2
    #
    b = np.arcsin(sin_b)
    l_NGP_l = np.arctan2(cos_b_sin_l_NGP_l, cos_b_cos_l_NGP_l)
    l = l_NGP - l_NGP_l
    l = (l + np.pi) % (2 * np.pi) - np.pi
    return l, b


def galactic2equatorial(l, b):
    sin_delta_1 = np.sin(delta_NGP) * np.sin(b)
    sin_delta_2 = np.cos(delta_NGP) * np.cos(b) * np.cos(l_NGP - l)
    sin_delta = sin_delta_1 + sin_delta_2
    cos_delta_sin_a_a_NGP = np.cos(b) * np.sin(l_NGP - l)
    cos_delta_cos_a_a_NGP_1 = np.cos(delta_NGP) * np.sin(b)
    cos_delta_cos_a_a_NGP_2 = np.sin(delta_NGP) * np.cos(b) * np.cos(l_NGP - l)
    cos_delta_cos_a_a_NGP = cos_delta_cos_a_a_NGP_1 - cos_delta_cos_a_a_NGP_2
    #
    delta = np.arcsin(sin_delta)
    a_a_NGP = np.arctan2(cos_delta_sin_a_a_NGP, cos_delta_cos_a_a_NGP)
    alpha = a_a_NGP + alpha_NGP
    alpha = (alpha + np.pi) % (2 * np.pi) - np.pi
    return alpha, delta
