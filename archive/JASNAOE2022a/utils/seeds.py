import numpy as np


def gen_seeds(master_seed, n):
    np.random.seed(seed=master_seed)
    seeds = np.random.randint(0, 2**31, n)
    return seeds
