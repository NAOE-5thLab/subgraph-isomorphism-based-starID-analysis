from matplotlib.pyplot import get
import numpy as np


def rodrigues_rotation_matrix(theta, n):
    n1 = n[0]
    n2 = n[1]
    n3 = n[2]
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    #
    R = np.array([
        [n1**2*(1-cos_theta) + cos_theta, n1*n2*(1-cos_theta) -
         n3*sin_theta, n1*n3*(1-cos_theta) + n2*sin_theta],
        [n1*n2*(1-cos_theta) + n3*sin_theta, n2**2*(1-cos_theta) +
         cos_theta, n2*n3*(1-cos_theta) - n1*sin_theta],
        [n1*n3*(1-cos_theta) - n2*sin_theta, n2*n3*(1-cos_theta) +
         n1*sin_theta, n3**2*(1-cos_theta) + cos_theta]
    ])
    return R


def rotation_on_earth(vec, alpha, delta, A, get_basis=False):
    e_x = np.array([1, 0, 0])
    e_y = np.array([0, 1, 0])
    e_z = np.array([0, 0, 1])
    #
    e_x_alpha = rodrigues_rotation_matrix(alpha, e_z) @ e_x
    e_y_alpha = rodrigues_rotation_matrix(alpha, e_z) @ e_y
    e_z_alpha = rodrigues_rotation_matrix(alpha, e_z) @ e_z
    vec_alpha = rodrigues_rotation_matrix(alpha, e_z) @ vec
    #
    e_x_alpha_delta = rodrigues_rotation_matrix(delta, -e_y_alpha) @ e_x_alpha
    e_y_alpha_delta = rodrigues_rotation_matrix(delta, -e_y_alpha) @ e_y_alpha
    e_z_alpha_delta = rodrigues_rotation_matrix(delta, -e_y_alpha) @ e_z_alpha
    vec_alpha_delta = rodrigues_rotation_matrix(delta, -e_y_alpha) @ vec_alpha
    #
    e_x_alpha_delta_A = rodrigues_rotation_matrix(
        A, -e_x_alpha_delta) @ e_x_alpha_delta
    e_y_alpha_delta_A = rodrigues_rotation_matrix(
        A, -e_x_alpha_delta) @ e_y_alpha_delta
    e_z_alpha_delta_A = rodrigues_rotation_matrix(
        A, -e_x_alpha_delta) @ e_z_alpha_delta
    vec_alpha_delta_A = rodrigues_rotation_matrix(
        A, -e_x_alpha_delta) @ vec_alpha_delta
    if get_basis:
        info = {
            'e_x': e_x,
            'e_y': e_y,
            'e_z': e_z,
            'e_x_alpha': e_x_alpha,
            'e_y_alpha': e_y_alpha,
            'e_z_alpha': e_z_alpha,
            'e_x_alpha_delta': e_x_alpha_delta,
            'e_y_alpha_delta': e_y_alpha_delta,
            'e_z_alpha_delta': e_z_alpha_delta,
            'e_x_alpha_delta_A': e_x_alpha_delta_A,
            'e_y_alpha_delta_A': e_y_alpha_delta_A,
            'e_z_alpha_delta_A': e_z_alpha_delta_A,
        }
        return vec_alpha_delta_A, info
    else:
        return vec_alpha_delta_A
