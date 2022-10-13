import numpy as np


def horizontal2vec(A, h):
    x = np.sin(h)
    y = np.sin(A)*np.cos(h)
    z = np.cos(A)*np.cos(h)
    return np.array([x, y, z])