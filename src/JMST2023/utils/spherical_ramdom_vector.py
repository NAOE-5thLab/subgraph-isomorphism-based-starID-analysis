import numpy as np
from utils.seeds import gen_seeds
from utils.soler import newton_method, bisection_method


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
    b = np.concatenate(
        [x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]],
        axis=1)
    if n == 1:
        b = b[0]
    return b


def limit_uniform_spherical_vector(sigma=10*np.pi/180, n=1, seed=100):
    seeds = gen_seeds(seed, 2)
    np.random.seed(seed=seeds[0])
    u = np.random.rand(n)
    np.random.seed(seed=seeds[1])
    v = np.random.rand(n)
    #
    z = 1 - (1 - np.cos(sigma))*u
    x = np.sqrt(1-z**2)*np.cos(2*np.pi*v)
    y = np.sqrt(1-z**2)*np.sin(2*np.pi*v)
    b = np.concatenate(
        [x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]],
        axis=1)
    if n == 1:
        b = b[0]
    return b


def normal_spherical_approximated_vector(mu, sigma, n=1, seed=100):
    eps = 1.0e-8
    if np.linalg.norm(mu) - 1.0 > eps:
        print('!!!caution!!! : The length of input vector is not 1.')
        mu = mu/np.linalg.norm(mu)
    #
    n_z_temp = np.array([-mu[2]*mu[0], -mu[2]*mu[1], 1 - mu[2]**2])
    n_z = n_z_temp/np.linalg.norm(n_z_temp)
    n_y = np.cross(n_z, mu)
    #
    seeds = gen_seeds(seed, 2)
    np.random.seed(seed=seeds[0])
    k_y = np.random.randn(n)*sigma
    np.random.seed(seed=seeds[1])
    k_z = np.random.randn(n)*sigma
    vectors_hat_temp = mu + k_y*n_y + k_z*n_z
    if n > 1:
        vectors_hat = vectors_hat_temp/np.linalg.norm(vectors_hat_temp, axis=1)
    else:
        vectors_hat = vectors_hat_temp/np.linalg.norm(vectors_hat_temp)
    return vectors_hat


def normal_spherical_vector(mu, kappa, n=1, seed=100):
    rand = VonMisesFisherDistribution(mu, kappa)
    return rand.sampling(n, seed=seed)


class VonMisesFisherDistribution:
    def __init__(self, mu, kappa):
        self.mu = mu
        self.kappa = kappa
        self.inv_kappa = 1/self.kappa
        self.exp_m2kappa = np.exp(-2 * self.kappa)
        #
        inv_norm = 1/np.linalg.norm(self.mu)
        mx = self.mu[0]*inv_norm
        my = self.mu[1]*inv_norm
        mz = self.mu[2]*inv_norm
        #
        angle = np.arccos(mz)
        c = np.cos(angle * 0.5)
        s = np.sin(angle * 0.5)
        norm = np.sqrt(mx * mx + my * my)
        #
        if norm > 0:
            qr = c
            qi = s * -my / norm
            qj = s * mx / norm
        else:
            qr = 1
            qi = 0
            qj = 0
        self.qri = qr * qi
        self.qij = qi * qj
        self.qrj = qr * qj
        self.qxx = qr * qr + qi * qi - qj * qj
        self.qyy = qr * qr - qi * qi + qj * qj
        self.qzz = qr * qr - qi * qi - qj * qj

    def sampling(self, n=1, seed=100):
        seeds = gen_seeds(seed, 2)
        np.random.seed(seed=seeds[0])
        r = np.random.rand(n)
        if self.kappa > 0:
            w = 1 + np.log(r + (1 - r) * self.exp_m2kappa) * self.inv_kappa
        else:
            w = 2 * r - 1
        s = np.sqrt(1 - w * w)
        np.random.seed(seed=seeds[1])
        theta = 2 * np.pi * np.random.rand(n)
        #
        x = s * np.cos(theta)
        y = s * np.sin(theta)
        z = w
        #
        xx = x * self.qxx + 2 * (y * self.qij + z * self.qrj)
        yy = y * self.qyy - 2 * (z * self.qri - x * self.qij)
        zz = z * self.qzz - 2 * (x * self.qrj - y * self.qri)
        sample = np.concatenate(
            [xx[:, np.newaxis], yy[:, np.newaxis], zz[:, np.newaxis]], axis=1)
        if n == 1:
            sample = sample[0]
        return sample


def estimate_kappa(theta, alpha, exp_flag=False):
    def f(kappa_temp):
        if exp_flag:
            kappa = np.exp(kappa_temp)
            return von_mises_fisher_cum_dist_func(theta, kappa) - alpha
        else:
            kappa = kappa_temp
            return von_mises_fisher_cum_dist_func(theta, kappa) - alpha
    if exp_flag:
        kappa_min = 1.0e-16
        kappa_max = 50
        est_kappa_temp = bisection_method(f, kappa_min, kappa_max)
        est_kappa = np.exp(est_kappa_temp)
        return est_kappa
    else:
        return newton_method(f, 1.0)


def von_mises_fisher_cum_dist_func(theta, kappa):
    numerator = np.exp(-2*kappa) - np.exp(kappa * (np.cos(theta) - 1.0))
    denominator = 1 - np.exp(-2*kappa)
    return 1 + numerator/denominator
