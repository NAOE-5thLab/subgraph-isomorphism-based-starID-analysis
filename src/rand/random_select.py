import random
from select import select
import numpy as np

from rand.spherical_uniform import uniform_spherical_vector


def random_select(n, index, seed=100):
    random.seed(seed)
    selected_index = random.sample(index, n)
    return selected_index

