import numpy as np

from utils.seeds import gen_seeds
from utils.quaternion import convert_z_axis


def uniform_spherical_vector(n=1, seed=100):
    seeds = gen_seeds(seed, 2)
    np.random.seed(seed=seeds[0])
    u = np.random.rand(n)
    np.random.seed(seed=seeds[1])
    v = np.random.rand(n)
    #
    z = -2*u + 1
    x = np.sqrt(1-z**2)*np.cos(2*np.pi*v)
    y = np.sqrt(1-z**2)*np.sin(2*np.pi*v)
    sample = np.concatenate([
        x[:, np.newaxis],
        y[:, np.newaxis],
        z[:, np.newaxis]],
        axis=1)
    if n == 1:
        sample = sample[0]
    return sample


def uniform_limited_spherical_vector(mu, theta_max, n=1, seed=10, theta_flag=False):
    seeds = gen_seeds(seed, 2)
    np.random.seed(seed=seeds[0])
    u = np.random.rand(n)
    np.random.seed(seed=seeds[1])
    v = np.random.rand(n)
    #
    z = 1 - (1 - np.cos(theta_max)) * u
    x = np.sqrt(1-z**2)*np.cos(2*np.pi*v)
    y = np.sqrt(1-z**2)*np.sin(2*np.pi*v)
    vec = np.concatenate([
        x[:, np.newaxis],
        y[:, np.newaxis],
        z[:, np.newaxis]],
        axis=1)
    sample = convert_z_axis(vec, mu)
    if n == 1:
        sample = sample[0]
    if theta_flag:
        return sample, np.arccos(z)
    return sample
