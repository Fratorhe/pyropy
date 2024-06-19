from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

COL_NAMES = ["T", "time", "DSC", "mass", "gas", "sens.", "segment"]


##Temp./�C;Time/min;DSC/(mW/mg);Mass/%;Gas Flow(protective)/(ml/min);Sensit./(uV/mW);Segment


@dataclass
class ExperimentReaderNetzsch:
    filename: str
    rate: float = 20  # K/min
    folder: str = "./"
    segment: int = 2
    savgol_win_length: int = 50
    savgol_poly_order: int = 2

    def __post_init__(self):
        # Read the file, skipping the first 36 lines
        df = pd.read_csv(
            f"{self.folder}/{self.filename}", skiprows=37, delimiter=";", encoding="latin1", names=COL_NAMES
        )

        self.rateKs = self.rate / 60
        # convert T to K, time to s
        df["T"] = df["T"] + 273.15
        df["time"] = df["time"] * 60

        # select segment
        df = df[df["segment"] == 2].reset_index(drop=True)
        df["time"] = df["time"] - df["time"].iloc[0]  # reset the time
        df["mass"] = df["mass"] - df["mass"].iloc[0]  # reset the mass
        # smoothen it for better derivative computation
        df["mass"] = savgol_filter(
            df["mass"], window_length=self.savgol_win_length, polyorder=self.savgol_poly_order
        )

        # compute the derivative
        # we do derivative wrt time and divide by heating rate applying chain rule
        # dm/dT = dm/dt * 1/beta
        # we do it this way because T can oscilate, while t always increases.
        df["dm"] = -self.compute_derivative(df["time"].values, df["mass"].values) / self.rateKs

        ## if need to apply the filter in the derivative use this:
        # Apply the Savitzky-Golay filter to smooth the derivative
        # dy_dx_smooth = savgol_filter(df['dm'], window_length=31, polyorder=2)

        # Add the smoothed derivative to the DataFrame
        # df['dm'] = dy_dx_smooth

        self.time = df["time"].values
        self.temperature = df["T"].values
        self.dRho = df["dm"].values
        self.Rho = df["mass"].values

    def plotit(self):
        fig, (tga, dtga) = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

        tga.plot(self.temperature, self.Rho)
        dtga.plot(self.temperature, self.dRho)

        tga.set_ylabel("TGA")

        dtga.set_xlabel("Time")
        dtga.set_ylabel("dTGA")

    def compute_derivative(self, x, y):
        # Compute the derivative using finite differences
        dy_dx = np.zeros_like(y)

        # Central differences for the interior points
        dy_dx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])

        # Forward difference for the first point
        dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])

        # Backward difference for the last point
        dy_dx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

        # Add the derivative to the DataFrame
        return dy_dx
