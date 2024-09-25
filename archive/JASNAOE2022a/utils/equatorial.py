import numpy as np


def vec2equatorial(vec):
    if vec.shape == (3,):
        x = vec[0]
        y = vec[1]
        z = vec[2]
    else:
        x = vec[:, 0]
        y = vec[:, 1]
        z = vec[:, 2]
    # 
    alpha = np.arctan(y/x)
    # 
    if vec.shape == (3,):
        if (x<0)*(y>=0):
            alpha += np.pi
        elif (x<0)*(y<0):
            alpha -= np.pi
    else:
        alpha[(x<0)*(y>=0)] += np.pi
        alpha[(x<0)*(y<0)] -= np.pi
    delta = np.pi/2 - np.arccos(z)
    return alpha, delta


def equatorial2vec(alpha, delta):
    x = np.cos(alpha) * np.cos(delta)
    y = np.sin(alpha) * np.cos(delta)
    z = np.sin(delta)
    if alpha.shape == ():
        return np.array([x, y, z])
    else:
        return np.concatenate(
            [x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]], axis=1)
