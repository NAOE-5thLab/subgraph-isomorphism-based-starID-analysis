import sys
import numpy as np
import pandas as pd
import tqdm


EPS = sys.float_info.epsilon


def get_j_by_binary_search(vector, value):
    left_index = 0
    right_index = len(vector) - 1
    #
    center_index = int((left_index + right_index)/2 + 0.5)
    while right_index - left_index > 1:
        if value <= vector[center_index]:
            right_index = center_index
        else:
            left_index = center_index
        center_index = int((left_index + right_index)/2 + 0.5)
    return left_index


def check_fname(fname, extension='.csv'):
    if len(fname) < len(extension):
        fname += extension
    else:
        if fname[len(extension):] != extension:
            fname += extension
    return fname


class IndexSearch:
    def __init__(self, y_vector=None):
        self.y_vector = np.array(y_vector)
        if y_vector is None:
            print('!!! Warnning : You should load data !!!')
        else:
            self.i_vector = np.argsort(self.y_vector)
            self.s_vector = self.y_vector[self.i_vector]

    def get_interval_index(self, y_lower_bound, y_upper_bound):
        lower_erase = y_lower_bound < self.s_vector
        upper_erase = self.s_vector < y_upper_bound
        clip = lower_erase*upper_erase
        interval_index = self.i_vector[clip]
        return interval_index

    def save_vector(self, dir, fname):
        df = pd.DataFrame(
            [self.y_vector, self.i_vector, self.s_vector],
            columns=['y', 'i', 's'])
        fname = check_fname(fname, extension='.csv')
        df.to_csv(dir + fname)

    def load_vector(self, dir, fname):
        fname = check_fname(fname, extension='.csv')
        df = pd.read_csv(dir + fname)
        self.y_vector = df['y'].to_numpy()
        self.i_vector = df['i'].to_numpy()
        self.s_vector = df['s'].to_numpy()


class KVector(IndexSearch):
    def __init__(self, y_vector=None):
        super().__init__(y_vector)
        self.k_vector = None

    def create_kvector(self):
        self.N = len(self.y_vector)
        self.y_max = np.max(self.y_vector)
        self.y_min = np.min(self.y_vector)
        self.xi = EPS*max(abs(self.y_max), abs(self.y_min))
        self.m = (self.y_max - self.y_min + 2*self.xi)/(self.N - 1)
        self.q = self.y_min - self.m - self.xi
        #
        if self.k_vector is None:
            print('Creating k-vector')
            k_vector = []
            for i in tqdm.tqdm(range(self.N)):
                if i == 0:
                    k_vector.append(0)
                elif i == self.N - 1:
                    k_vector.append(self.N)
                else:
                    x = i + 1
                    z = self.m*x + self.q
                    j = get_j_by_binary_search(self.s_vector, z)
                    k_vector.append(j)
            self.k_vector = np.array(k_vector)

    def get_interval_index(self, y_lower_bound, y_upper_bound):
        jb = int((y_lower_bound - self.q)/self.m)
        jt = int((y_upper_bound - self.q)/self.m + 1.0)
        #
        k_index_start = self.k_vector[jb-1] + 1
        k_index_end = self.k_vector[jt-1]
        #
        interval_index = self.i_vector[k_index_start-1:k_index_end]
        return interval_index

    def save_vector(self, path):
        df = pd.DataFrame(
            [self.y_vector, self.i_vector, self.s_vector, self.k_vector],
            columns=['y', 'i', 's', 'k'])
        df.to_csv(path + '.csv')

    def load_vector(self, path):
        df = pd.read_csv(path + '.csv')
        self.y_vector = df['y'].to_numpy()
        self.i_vector = df['i'].to_numpy()
        self.s_vector = df['s'].to_numpy()
        self.k_vector = df['k'].to_numpy()


# def main():
#     pass


# if __name__ == '__main__':
#     main()
#     print('tests were completed')
