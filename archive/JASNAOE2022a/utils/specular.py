import numpy as np


def specular_sign(vec1, vec2, vec3):
    cross = np.cross(vec2, vec3)
    inner = np.inner(vec1, cross)
    return np.sign(inner)
