import numpy as np

# from numba import jit


# @jit
# def vec2equatorial(vec):
#     e_x = np.array([1, 0, 0])
#     e_y = np.array([0, 1, 0])
#     e_z = np.array([0, 0, 1])
#     #
#     vec_dot_ez = np.inner(vec, e_z)
#     vec_abs = np.linalg.norm(vec)
#     ez_abs = np.linalg.norm(e_z)
#     delta = np.pi / 2 - np.arccos(vec_dot_ez / (vec_abs * ez_abs))
#     #
#     orthographic_ez_s1 = vec_dot_ez / ez_abs**2 * e_z
#     orthographic_ez_s1_prime = vec - orthographic_ez_s1
#     o_dot_ex = np.inner(orthographic_ez_s1_prime, e_x)
#     o_abs = np.linalg.norm(orthographic_ez_s1_prime)
#     ex_abs = np.linalg.norm(e_x)
#     alpha = np.arccos(o_dot_ex / (o_abs * ex_abs))
#     return alpha, delta


def vec2equatorial(vec):
    vec = vec / np.linalg.norm(vec)
    delta = np.arcsin(vec[2])
    alpha = np.arctan2(vec[1], vec[0])
    alpha = (alpha + np.pi) % (2 * np.pi) - np.pi
    return alpha, delta


# @jit
def equatorial2vec(alpha, delta):
    if isinstance(alpha, (int, float)):
        alpha = np.array([alpha])
        delta = np.array([delta])
        rets = np.empty((1, 3))
    else:
        rets = np.empty((len(alpha), 3))

    rets[:, 0] = np.cos(alpha) * np.cos(delta)
    rets[:, 1] = np.sin(alpha) * np.cos(delta)
    rets[:, 2] = np.sin(delta)

    if rets.shape[0] == 1:
        return rets[0]
    else:
        return rets
