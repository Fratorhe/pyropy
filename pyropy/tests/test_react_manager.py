import os

import pytest

from pyropy import ReactManager

file_path = os.path.dirname(__file__)


@pytest.fixture
def reaction_manager() -> ReactManager:
    return ReactManager(filename="data_parallel.json", folder=f"{file_path}/")


def test_number_reactions(reaction_manager: ReactManager) -> None:
    reaction_manager.react_reader()
    assert reaction_manager.n_reactions == 2


def test_number_solids(reaction_manager: ReactManager) -> None:
    reaction_manager.react_reader()
    assert reaction_manager.n_solids == 1


def test_number_solids(reaction_manager: ReactManager) -> None:
    reaction_manager.react_reader()
    assert len(reaction_manager.gas_product) == 2
