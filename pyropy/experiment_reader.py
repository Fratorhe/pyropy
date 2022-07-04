from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd


class ExperimentReader(Protocol):
    """
    This protocol specifies how a ExperimentReader should look like.
    It should take as inputs
    """

    filename: str
    folder: str
    time: np.array
    temperature: np.array
    dRho: np.array
    Rho: np.array


@dataclass
class ExperimentReaderCSV(ExperimentReader):
    """ """

    filename: str
    folder: str = "./"

    def __post_init__(self):
        file = pd.read_csv(self.folder + "/" + self.filename)
        self.time = file["time"].values
        self.temperature = file["temperature"].values
        self.dRho = file["dRho"].values
        self.Rho = file["rho"].values
