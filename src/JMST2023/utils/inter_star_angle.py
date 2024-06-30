import numpy as np
# from numba import jit


# @jit
def inter_star_angle_RADE(alpha1, delta1, alpha2, delta2):
    """inter star angle

    Equation (using spherical trigonometry)
    cos(d) = cos(pi/2 - delta1) * cos(pi/2 - delta2)
        + sin(pi/2 - delta1) * sin(pi/2 - delta2) * cos(alpha1 - delta2)

    Args:
        alpha1 (float): Right Ascension of star1
        delta1 (float): Declination of star1
        alpha2 (float): Right Ascension of star2
        delta2 (float): Declination of star2

    Returns:
        float: inter star angle
    """
    cos_inter_angle = np.sin(delta1) * np.sin(delta2) + np.cos(delta1) * \
        np.cos(delta2) * np.cos(alpha1 - alpha2)
    inter_angle = np.arccos(cos_inter_angle)
    return inter_angle


# @jit
def inter_star_angle_vec(star1_vec, star2_vec):
    inner = np.inner(star1_vec, star2_vec)
    star1_len = np.linalg.norm(star1_vec, axis=star1_vec.ndim-1)
    star2_len = np.linalg.norm(star2_vec, axis=star2_vec.ndim-1)
    cos = inner/(star1_len*star2_len)
    theta = np.arccos(cos)
    return theta
