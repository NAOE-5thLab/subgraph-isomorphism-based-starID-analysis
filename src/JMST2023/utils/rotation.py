from matplotlib.pyplot import get
import numpy as np


def rodrigues_rotation_matrix(theta, n):
    if n.ndim == 1:
        n1 = n[0]
        n2 = n[1]
        n3 = n[2]
    elif n.ndim == 2:
        n1 = n[:, 0]
        n2 = n[:, 1]
        n3 = n[:, 2]
    else:
        print('!!! Warning : input dimension error !!!')
        return 0
    #
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    #
    r11 = n1**2*(1-cos_theta) + cos_theta
    r12 = n1*n2*(1-cos_theta) - n3*sin_theta
    r13 = n1*n3*(1-cos_theta) + n2*sin_theta
    r21 = n1*n2*(1-cos_theta) + n3*sin_theta
    r22 = n2**2*(1-cos_theta) + cos_theta
    r23 = n2*n3*(1-cos_theta) - n1*sin_theta
    r31 = n1*n3*(1-cos_theta) - n2*sin_theta
    r32 = n2*n3*(1-cos_theta) + n1*sin_theta
    r33 = n3**2*(1-cos_theta) + cos_theta
    #
    if r11.shape == ():
        R = np.array([
            [r11, r12, r13],
            [r21, r22, r23],
            [r31, r32, r33]
        ])
    else:
        r1 = np.concatenate(
            [r11[:, np.newaxis, np.newaxis], r12[:, np.newaxis, np.newaxis], r13[:, np.newaxis, np.newaxis]], axis=2)
        r2 = np.concatenate(
            [r21[:, np.newaxis, np.newaxis], r22[:, np.newaxis, np.newaxis], r23[:, np.newaxis, np.newaxis]], axis=2)
        r3 = np.concatenate(
            [r31[:, np.newaxis, np.newaxis], r32[:, np.newaxis, np.newaxis], r33[:, np.newaxis, np.newaxis]], axis=2)
        R = np.concatenate(
            [r1, r2, r3], axis=1)
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
