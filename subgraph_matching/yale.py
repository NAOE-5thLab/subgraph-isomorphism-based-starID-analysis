import os
import numpy as np
import pandas as pd


class YaleStarCatalog:
    """
    Class to load and process the Yale Bright Star Catalog.
    """

    # File Path
    data_path = f"{os.path.dirname(__file__)}/catalogue_gz/bsc5.dat"

    # Parameters determined by './catalogue_gz/readme.txt'
    columns = (
        ["HR", "Name", "DM", "HD", "SAO"]
        + ["FK5", "IRflag", "r_IRflag*"]
        + ["Multiple*", "ADS", "ADScomp"]
        + ["VarID", "RAh1900 [h]", "RAm1900 [min]"]
        + ["RAs1900 [s]", "DE-1900", "DEd1900 [deg]"]
        + ["DEm1900 [arcmin]", "DEs1900 [arcsec]"]
        + ["RAh [h]", "RAm [min]", "RAs [s]", "DE-"]
        + ["DEd [deg]", "DEm [arcmin]", "DEs [arcsec]"]
        + ["GLON [deg]", "GLAT [deg]", "Vmag [mag]"]
        + ["n_Vmag*", "u_Vmag", "B-V [mag]", "u_B-V"]
        + ["U-B [mag]", "u_U-B", "R-I [mag]", "n_R-I"]
        + ["SpType", "n_SpType", "pmRA [arcsec/yr]*"]
        + ["pmDE [arcsec/yr]", "n_Parallax"]
        + ["Parallax [arcsec]", "RadVel [km/s]"]
        + ["n_RadVel*", "l_RotVel", "RotVel [km/s]"]
        + ["u_RotVel", "Dmag [mag]", "Sep [arcsec]"]
        + ["MultID", "MultCnt", "NoteFlag"]
    )
    dtype_label = (
        ["I4", "A10", "A11", "I6", "I6", "I4", "A1"]
        + ["A1", "A1", "A5", "A2", "A9", "I2", "I2"]
        + ["F4.1", "A1", "I2", "I2", "I2", "I2", "I2"]
        + ["F4.1", "A1", "I2", "I2", "I2", "F6.2", "F6.2"]
        + ["F5.2", "A1", "A1", "F5.2", "A1", "F5.2", "A1"]
        + ["F5.2", "A1", "A20", "A1", "F6.3", "F6.3", "A1"]
        + ["F5.3", "I4", "A4", "A2", "I3", "A1", "F4.1"]
        + ["F6.1", "A4", "I2", "A1"]
    )
    split_id = (
        [0, 4, 14, 25, 31, 37, 41, 42]
        + [43, 44, 49, 51, 60, 62, 64, 68]
        + [69, 71, 73, 75, 77, 79, 83, 84]
        + [86, 88, 90, 96, 102, 107, 108]
        + [109, 114, 115, 120, 121, 126]
        + [127, 147, 148, 154, 160, 161]
        + [166, 170, 174, 176, 179, 180]
        + [184, 190, 194, 196, 197]
    )
    split_len = len(split_id)

    def __init__(self, log_dir="./"):
        """
        Initialize the YaleStarCatalog class.

        Args:
            log_dir (str): Path to the log directory.
        """
        self.log_dir = log_dir
        path = f"{self.log_dir}/temp.csv"
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        self.load()
        self.reshape()
        self.save_df()

    def get_info(self):
        """
        Print information about the catalog.
        """
        print("----- Yale Bright Star Catalog -----")
        Vmag = self.get_Vmag()
        print(f"the number of stars : {len(self.get_HR())}")
        print(f"the range of magnitude : [{min(Vmag)}, {max(Vmag)}]")

    def load(self):
        """
        Load the data file and convert it to a DataFrame.
        """
        with open(self.data_path, "r") as f:
            row_list = f.read().split("\n")

        row_data_list = []
        for row in row_list:
            row_data = []
            for i in range(self.split_len - 1):
                data = row[self.split_id[i] : self.split_id[i + 1]].strip()
                if self.dtype_label[i][0] == "I":
                    data = int(data) if data != "" else None
                elif self.dtype_label[i][0] == "F":
                    data = float(data) if data != "" else None
                row_data.append(data)
            row_data_list.append(row_data)
        self.row_df = pd.DataFrame(row_data_list, columns=self.columns)

    def reshape(self):
        """
        Reshape the DataFrame by removing missing values and converting units.
        """
        self.reshaped_df = self.row_df.copy()
        # Erase missing values
        self.reshaped_df = self.reshaped_df.dropna(subset=["Vmag [mag]"])
        # Convert units
        self.convert_units()
        # Convert HR to integer
        self.reshaped_df["HR"] = self.reshaped_df["HR"].astype("int")

    def convert_units(self):
        """
        Convert units for Right Ascension and Declination.
        """
        # Right Ascension in [0, 2pi]
        self.reshaped_df["RA [h]"] = (
            self.row_df["RAh [h]"]
            + self.row_df["RAm [min]"] / 60.0
            + self.row_df["RAs [s]"] / 60.0**2
        )
        self.reshaped_df["RA [rad]"] = self.reshaped_df["RA [h]"] * np.pi / 12

        # Declination in [-pi, pi]
        self.reshaped_df["DE [deg]"] = (
            self.row_df["DEd [deg]"]
            + self.row_df["DEm [arcmin]"] / 60.0
            + self.row_df["DEs [arcsec]"] / 60.0**2
        )
        minus_index = self.reshaped_df["DE-"] == "-"
        self.reshaped_df.loc[minus_index, "DE [deg]"] = -self.reshaped_df["DE [deg]"]
        self.reshaped_df["DE [rad]"] = self.reshaped_df["DE [deg]"] * np.pi / 180

    def save_df(self):
        """
        Save the reshaped DataFrame to a CSV file.
        """
        self.reshaped_df[["HR", "RA [rad]", "DE [rad]", "Vmag [mag]", "VarID"]].to_csv(
            f"{self.log_dir}/yale_catalog.csv"
        )

    def get_df(self):
        """
        Get the reshaped DataFrame with selected columns.

        Returns:
            DataFrame: The reshaped DataFrame.
        """
        return self.reshaped_df[["HR", "RA [rad]", "DE [rad]", "Vmag [mag]"]]

    def get_full_df(self):
        """
        Get the reshaped DataFrame.

        Returns:
            DataFrame: The reshaped DataFrame.
        """
        return self.reshaped_df

    def get_D(self):
        """
        Get the reshaped DataFrame as a numpy array.

        Returns:
            numpy.ndarray: The reshaped DataFrame as a numpy array.
        """
        return self.get_df().to_numpy()

    def get_HR(self):
        """
        Get the HR column as a numpy array.

        Returns:
            numpy.ndarray: The HR column.
        """
        return self.reshaped_df["HR"].to_numpy()

    def get_RA(self):
        """
        Get the Right Ascension column as a numpy array.

        Returns:
            numpy.ndarray: The Right Ascension column.
        """
        return self.reshaped_df["RA [rad]"].to_numpy()

    def get_DE(self):
        """
        Get the Declination column as a numpy array.

        Returns:
            numpy.ndarray: The Declination column.
        """
        return self.reshaped_df["DE [rad]"].to_numpy()

    def get_Vmag(self):
        """
        Get the Vmag column as a numpy array.

        Returns:
            numpy.ndarray: The Vmag column.
        """
        return self.reshaped_df["Vmag [mag]"].to_numpy()

    def get_VarID(self):
        """
        Get the VarID column as a numpy array.

        Returns:
            numpy.ndarray: The VarID column.
        """
        return self.reshaped_df["VarID"].to_numpy()


def main():
    db = YaleStarCatalog()
    print(type(db.get_HR()))


if __name__ == "__main__":
    main()
    print("task completed")
