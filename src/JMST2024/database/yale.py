import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(__file__)


class YaleStarCatalog:
    # File Path
    source_data_path = '/gz/'
    data_fname = 'bsc5.dat'
    note_fname = 'bsc5.notes'
    # Parameter determined by 'gz/readme.txt'
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

    def __init__(self, log_dir='./'):
        self.log_dir = log_dir
        print('----- Yale Bright Star Catalog -----')
        self.load()
        self.reshape()
        self.save_df()
        #
        Vmag = self.get_Vmag()
        print(f"the number of stars : {len(self.get_HR())}")
        print(f"the range of magnitude : [{min(Vmag)}, {max(Vmag)}]")

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
        self.row_df = pd.DataFrame(row_data_list, columns=self.columns)

    def reshape(self):
        self.reshaped_df = self.row_df.copy()
        # erase missing value
        self.reshaped_df = self.reshaped_df.dropna(subset=['Vmag [mag]'])
        # convert unit
        # Right Ascension in [0, 2pi]
        self.reshaped_df['RA [h]'] = self.row_df['RAh [h]'] + \
            self.row_df['RAm [min]']/60.0 + \
            self.row_df['RAs [s]']/60.0**2
        self.reshaped_df['RA [rad]'] = self.reshaped_df['RA [h]']*np.pi/12
        # Declination in [-pi, pi]
        self.reshaped_df['DE [deg]'] = self.row_df['DEd [deg]'] + \
            self.row_df['DEm [arcmin]']/60.0 + \
            self.row_df['DEs [arcsec]']/60.0**2
        minus_index = list(self.reshaped_df['DE-'] == '-')
        self.reshaped_df.loc[minus_index, 'DE [deg]'] = - \
            self.reshaped_df['DE [deg]']
        self.reshaped_df['DE [rad]'] = self.reshaped_df['DE [deg]']*np.pi/180
        # cnvert HR to integer
        self.reshaped_df['HR'] = self.reshaped_df['HR'].astype('int')

    def save_df(self):
        self.reshaped_df[['HR', 'RA [rad]', 'DE [rad]', 'Vmag [mag]', 'VarID']].to_csv(
            f'{self.log_dir}/yale_catalog.csv')

    def get_df(self):
        return self.reshaped_df[['HR', 'RA [rad]', 'DE [rad]', 'Vmag [mag]']]

    def get_HR(self):
        return self.reshaped_df['HR'].to_numpy()

    def get_RA(self):
        return self.reshaped_df['RA [rad]'].to_numpy()

    def get_DE(self):
        return self.reshaped_df['DE [rad]'].to_numpy()

    def get_Vmag(self):
        return self.reshaped_df['Vmag [mag]'].to_numpy()

    def get_VarID(self):
        return self.reshaped_df['VarID'].to_numpy()


def main():
    db = YaleStarCatalog()
    print(type(db.get_HR()))


if __name__ == "__main__":
    main()
    print('task completed')
