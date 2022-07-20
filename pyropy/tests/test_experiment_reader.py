import os

import pytest

from pyropy.experiment_reader import ExperimentReaderCSV

file_path = os.path.dirname(__file__)


def test_experiment_reader():
    """
    Test to check the file reader.
    In this case, we check that the number of rows is 8 as in the file data.csv
    :return:
    """
    experiment = ExperimentReaderCSV(filename="data.csv", folder=f"{file_path}/")
    assert len(experiment.temperature) == 8


def test_experiment_reader_file():
    """
    Test to check the file reader.
    In this case, we check that if the file does not exist, we should get an exception of FileNotFoundError
    """
    with pytest.raises(FileNotFoundError):
        experiment = ExperimentReaderCSV(filename="asdf")
