import os

import numpy as np
import pytest

from pyropy import ReactManager
from pyropy.pyrolysis import PyrolysisCompetitive


file_path = (os.path.dirname(__file__))

@pytest.fixture
def reactions() -> ReactManager:
    return ReactManager(filename="data_competing.json", folder=f'{file_path}/')


def test_competitive_scheme(reactions: ReactManager):
    reactions.react_reader()
    reactions.param_reader()
    test = PyrolysisCompetitive(temp_0=373, temp_end=400, beta=1, n_points=15, reaction_scheme_obj=reactions)
    # Numerical solution using 200 points
    test.solve_system()
    print(test.rho_solid)
    solution_known = np.array(
        [
            100.0,
            99.96653943,
            99.93004845,
            99.89031579,
            99.84710509,
            99.80017088,
            99.74925335,
            99.69408498,
            99.63438634,
            99.56986877,
            99.50023475,
            99.42517921,
            99.34439012,
            99.25755113,
            99.16434227,
        ]
    )

    np.testing.assert_allclose(test.rho_solid, solution_known)
