import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(__file__)


class Catalog:
    dataset_header = ['ID', 'RA [rad]', 'DE [rad]', 'Vmag [mag]']


class YaleStarCatalog(Catalog):
    # File Path
    source_data_path = '/db_yale/'
    data_fname = 'bsc5.dat'
    note_fname = 'bsc5.notes'
    # Parameter determined by 'db_yale/readme.txt'
    columns = [
        'HR', 'Name', 'DM', 'HD', 'SAO',
        'FK5', 'IRflag', 'r_IRflag*',
        'Multiple*', 'ADS', 'ADScomp',
        'VarID', 'RAh1900 [h]', 'RAm1900 [min]',
        'RAs1900 [s]', 'DE-1900', 'DEd1900 [deg]',
        'DEm1900 [arcmin]', 'DEs1900 [arcsec]',
        'RAh [h]', 'RAm [min]', 'RAs [s]', 'DE-',
        'DEd [deg]', 'DEm [arcmin]', 'DEs [arcsec]',
        'GLON [deg]', 'GLAT [deg]', 'Vmag [mag]',
        'n_Vmag*', 'u_Vmag', 'B-V [mag]', 'u_B-V',
        'U-B [mag]', 'u_U-B', 'R-I [mag]', 'n_R-I',
        'SpType', 'n_SpType', 'pmRA [arcsec/yr]*',
        'pmDE [arcsec/yr]', 'n_Parallax',
        'Parallax [arcsec]', 'RadVel [km/s]',
        'n_RadVel*', 'l_RotVel', 'RotVel [km/s]',
        'u_RotVel', 'Dmag [mag]', 'Sep [arcsec]',
        'MultID', 'MultCnt', 'NoteFlag']
    dtype_label = [
        'I4', 'A10', 'A11', 'I6', 'I6', 'I4', 'A1',
        'A1', 'A1', 'A5', 'A2', 'A9', 'I2', 'I2',
        'F4.1', 'A1', 'I2', 'I2', 'I2', 'I2', 'I2',
        'F4.1', 'A1', 'I2', 'I2', 'I2', 'F6.2', 'F6.2',
        'F5.2', 'A1', 'A1', 'F5.2', 'A1', 'F5.2', 'A1',
        'F5.2', 'A1', 'A20', 'A1', 'F6.3', 'F6.3', 'A1',
        'F5.3', 'I4', 'A4', 'A2', 'I3', 'A1', 'F4.1',
        'F6.1', 'A4', 'I2', 'A1'
    ]
    split_id = [
        0, 4, 14, 25, 31, 37, 41, 42,
        43, 44, 49, 51, 60, 62, 64, 68,
        69, 71, 73, 75, 77, 79, 83, 84,
        86, 88, 90, 96, 102, 107, 108,
        109, 114, 115, 120, 121, 126,
        127, 147, 148, 154, 160, 161,
        166, 170, 174, 176, 179, 180,
        184, 190, 194, 196, 197]
    split_len = len(split_id)

    def __init__(self):
        row_df = self.load()
        self.dataset = self.reshape(row_df)
        self.save_dataset(self.dataset)

    def load(self):
        with open(BASE_DIR + self.source_data_path + self.data_fname, 'r') as f:
            row_list = f.read().split("\n")
        #
        row_data_list = []
        for row in row_list:
            row_data = []
            for i in range(self.split_len-1):
                data = row[self.split_id[i]:self.split_id[i+1]].strip()
                if self.dtype_label[i][0] == 'I':
                    if data != '':
                        data = int(data)
                    else:
                        data = None
                elif self.dtype_label[i][0] == 'F':
                    if data != '':
                        data = float(data)
                    else:
                        data = None
                row_data.append(data)
            row_data_list.append(row_data)
        row_df = pd.DataFrame(row_data_list, columns=self.columns)
        return row_df

    def reshape(self, row_df):
        # create temporary box
        dataset_df = row_df.copy()
        # erase missing value
        dataset_df = dataset_df.dropna(subset=['Vmag [mag]'])

        ### convert unit ###
        # Right Ascension in [0, 2pi]
        dataset_df['RA [h]'] = row_df['RAh [h]'] + \
            row_df['RAm [min]']/60.0 + row_df['RAs [s]']/60.0**2
        dataset_df['RA [rad]'] = dataset_df['RA [h]']*np.pi/12
        # Declination in [-pi, pi]
        dataset_df['DE [deg]'] = row_df['DEd [deg]'] + \
            row_df['DEm [arcmin]']/60.0 + row_df['DEs [arcsec]']/60.0**2
        minus_index = list(dataset_df['DE-'] == '-')
        dataset_df.loc[minus_index, 'DE [deg]'] = -dataset_df['DE [deg]']
        dataset_df['DE [rad]'] = dataset_df['DE [deg]']*np.pi/180
        # cnvert HR to integer
        dataset_df['HR'] = dataset_df['HR'].astype('int')

        ### dataset ###
        dataset = dataset_df[
            ['HR', 'RA [rad]', 'DE [rad]', 'Vmag [mag]']].to_numpy()
        return dataset

    def get_I(self):
        return np.arange(len(self.dataset[:, 0]))

    def get_HR(self):
        return self.dataset[:, 0].astype('int64')

    def get_RA(self):
        return self.dataset[:, 1]

    def get_DE(self):
        return self.dataset[:, 2]

    def get_Vmag(self):
        return self.dataset[:, 3]

    def get_dataset(self):
        return self.dataset

    def info(self):
        print('----- Yale Bright Star Catalog -----')
        print(f"the number of stars : {len(self.get_ID())}")
        print(
            f"the range of magnitude : [{min(self.get_Vmag())}, {max(self.get_Vmag())}]")

    def save_dataset(self, dataset):
        np.savetxt(BASE_DIR + '/cache/dataset.csv', dataset, delimiter=',')


def main():
    db = YaleStarCatalog()
    print(type(db.get_HR()))


if __name__ == "__main__":
    main()
    print('task completed')
