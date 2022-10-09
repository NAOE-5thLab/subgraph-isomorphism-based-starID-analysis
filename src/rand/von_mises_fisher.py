import numpy as np

from utils.seeds import gen_seeds
from utils.solver import bisection_method
from utils.quaternion import convert_z_axis


def von_mises_fisher_3d_sampling(mu, kappa, n=1, seed=100, theta_flag=False):
    # uniform random number
    seeds = gen_seeds(seed, 2)
    np.random.seed(seed=seeds[0])
    u = np.random.rand(n)
    np.random.seed(seed=seeds[1])
    v = np.random.rand(n)
    # convert
    if kappa > 0:
        cos_theta = 1 + (1/kappa) * np.log(1 - (1 - np.exp(-2*kappa)) * u)
    else:
        cos_theta = 2 * u - 1
    sin_theta = np.sqrt(1 - cos_theta * cos_theta)
    phi = 2 * np.pi * v
    #
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta
    vec = np.concatenate([
        x[:, np.newaxis],
        y[:, np.newaxis],
        z[:, np.newaxis]],
        axis=1)
    sample = convert_z_axis(vec, mu)
    if n == 1:
        sample = sample[0]
        cos_theta = cos_theta[0]
    if theta_flag:
        return sample, np.arccos(cos_theta)
    return sample


def estimate_kappa_von_mises_fisher_3d(theta, alpha):
    def f(log_kappa):
        kappa = np.exp(log_kappa)
        numerator = np.exp(-2*kappa) - np.exp(kappa * (np.cos(theta) - 1.0))
        denominator = 1 - np.exp(-2*kappa)
        return 1 + numerator/denominator - alpha
    #
    log_kappa_min = -10
    log_kappa_max = 50
    est_log_kappa = bisection_method(f, log_kappa_min, log_kappa_max)
    est_kappa = np.exp(est_log_kappa)
    return est_kappa


def estimate_thata_von_mises_fisher_3d(kappa, alpha):
    def f(theta):
        numerator = np.exp(-2*kappa) - np.exp(kappa * (np.cos(theta) - 1.0))
        denominator = 1 - np.exp(-2*kappa)
        return 1 + numerator/denominator - alpha
    #
    log_kappa_min = 0
    log_kappa_max = np.pi
    est_theta = bisection_method(f, log_kappa_min, log_kappa_max)
    return est_theta
