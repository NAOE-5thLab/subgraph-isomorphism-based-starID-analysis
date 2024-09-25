import numpy as np
# from numba import jit


# @jit
def vec2equatorial(vec):
    e_x = np.array([1, 0, 0])
    e_y = np.array([0, 1, 0])
    e_z = np.array([0, 0, 1])
    #
    vec_dot_ez = np.inner(vec, e_z)
    vec_abs = np.linalg.norm(vec)
    ez_abs = np.linalg.norm(e_z)
    delta = np.pi/2 - np.arccos(vec_dot_ez/(vec_abs*ez_abs))
    #
    orthographic_ez_s1 = vec_dot_ez/ez_abs**2 * e_z
    orthographic_ez_s1_prime = vec - orthographic_ez_s1
    o_dot_ex = np.inner(orthographic_ez_s1_prime, e_x)
    o_abs = np.linalg.norm(orthographic_ez_s1_prime)
    ex_abs = np.linalg.norm(e_x)
    alpha = np.arccos(o_dot_ex/(o_abs*ex_abs))
    return alpha, delta


# @jit
def equatorial2vec(alpha, delta):
    x = np.cos(alpha) * np.cos(delta)
    y = np.sin(alpha) * np.cos(delta)
    z = np.sin(delta)
    if alpha.shape == ():
        return np.array([x, y, z])
    else:
        return np.concatenate(
            [x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]], axis=1)
