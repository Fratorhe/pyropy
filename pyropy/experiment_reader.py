from typing import Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd


class ExperimentReader(Protocol):
    """
    Protocol for experiment data readers.

    An implementing class should define attributes for file location and
    arrays containing experimental results.
    """

    def __init__(self, filename: str, folder: str) -> None:
        self.filename: str
        self.folder: str
        self.time: npt.ArrayLike
        self.temperature: npt.ArrayLike
        self.dRho: npt.ArrayLike
        self.Rho: npt.ArrayLike


class ExperimentReaderCSV:
    """
    CSV-based implementation of the ExperimentReader protocol.

    Reads experimental data from a CSV file containing columns:
    'time', 'temperature', 'dRho', and 'rho'.

    Attributes
    ----------
    filename : str
        Name of the CSV file.
    folder : str, optional
        Folder containing the file (default: './').
    time : np.ndarray
        Time data.
    temperature : np.ndarray
        Temperature data.
    dRho : np.ndarray
        Rate of change of density.
    Rho : np.ndarray
        Density data.
    """

    def __init__(self, filename: str, folder: str = "./") -> None:
        # Initialize attributes according to the protocol
        self.filename = filename
        self.folder = folder

        # Read the CSV file
        file = pd.read_csv(f"{self.folder}/{self.filename}")

        # Assign the data to the attributes defined in the protocol
        self.time = file["time"].to_numpy()
        self.temperature = file["temperature"].to_numpy()
        self.dRho = file["dRho"].to_numpy()
        self.Rho = file["rho"].to_numpy()
