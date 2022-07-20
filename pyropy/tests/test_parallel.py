import os

import numpy as np
import pytest

from pyropy import PyrolysisParallel, PyrolysisParallelAnalytical
from pyropy import ReactManager

file_path = (os.path.dirname(__file__))

@pytest.fixture
def reactions() -> ReactManager:
    return ReactManager(filename="data_parallel.json", folder=f'{file_path}/')


def test_parallel_scheme(reactions: ReactManager):
    T_0 = 373
    T_end = 2000
    b = 20
    reactions.react_reader()
    reactions.param_reader()
    test = PyrolysisParallel(temp_0=T_0, temp_end=T_end, beta=b, n_points=15, reaction_scheme_obj=reactions)
    # Numerical solution using 200 points
    test.solve_system()
    solution_known = np.array(
        [
            100.0,
            99.98303254,
            99.41627862,
            96.87519857,
            95.41962767,
            95.00241108,
            93.79111514,
            90.68643543,
            88.25904679,
            86.96281349,
            86.26421096,
            85.86255122,
            85.61696297,
            85.45879965,
            85.35238839,
        ]
    )

    np.testing.assert_allclose(test.rho_solid, solution_known)


def test_parallel_scheme_analytical(reactions: ReactManager):
    T_0 = 373
    T_end = 2000
    b = 20
    reactions.react_reader()
    reactions.param_reader()
    test = PyrolysisParallelAnalytical(temp_0=T_0, temp_end=T_end, beta=b, n_points=15, reaction_scheme_obj=reactions)
    # Numerical solution using 200 points
    test.solve_system()
    solution_known = np.array(
        [
            100.0,
            99.98303254,
            99.41627862,
            96.87519857,
            95.41962767,
            95.00241108,
            93.79111514,
            90.68643543,
            88.25904679,
            86.96281349,
            86.26421096,
            85.86255122,
            85.61696297,
            85.45879965,
            85.35238839,
        ]
    )

    np.testing.assert_allclose(test.rho_solid, solution_known, rtol=1e-05)

