from abc import ABC
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .reaction_reader_writer import ReactManager

R = 8.314


class Pyrolysis(Protocol):

    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False
    rho_solid: np.array([])
    drho_solid: np.array([])
    time: np.array([])
    temperature: np.array([])

    def generate_rate(self) -> list:
        ...

    def pyro_rates(self) -> list:
        ...

    def solve_system(self) -> None:
        ...

    def to_csv(self, filename: str) -> None:
        file_save = pd.DataFrame(data=None, columns=[], index=None)
        file_save['time'] = self.time
        file_save['temperature'] = self.temperature
        file_save['rho'] = self.rho_solid
        file_save['dRho'] = self.drho_solid
        file_save.to_csv(filename)


def compute_time(betaKs: float, temp_0: float = 373, temp_end: float = 2000, n_points: int = 500):
    """Compute temperature as a function of time based on the given heating rate and
    initial temperature temp_0. If no time vector is provided, a time vector of size n_points
    is computed from temp_0 to temp_end"""

    time = np.linspace(0, (temp_end - temp_0) / betaKs, n_points)
    return time

def compute_temperature(time: np.array, temp_0: float, betaKs: float):
    """Compute temperature as a function of time based on the given heating rate and
    initial temperature temp_0. If no time vector is provided, a time vector of size n_points
    is computed from temp_0 to temp_end"""

    temperature = time * betaKs + temp_0
    return temperature


@dataclass
class PyrolysisParallel(Pyrolysis):
    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False

    param_names: list = field(default_factory=list, init=False)
    dict_params: dict = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        # Convert beta in K/min to betaKs in K/sec
        self.betaKs = self.beta / 60
        # Compute temperature and time
        self.time = compute_time(temp_0=self.temp_0, temp_end=self.temp_end, betaKs=self.betaKs, n_points=self.n_points)
        self.temperature = compute_temperature(self.time, temp_0=self.temp_0, betaKs=self.betaKs)

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def generate_rate(self, temperature=float("inf")):
        k = []
        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            k.append(
                10 ** self.reaction_scheme_obj.dict_params["A"][idx]
                * np.exp(-self.reaction_scheme_obj.dict_params["E"][idx] / (R * temperature))
            )
        return k

    def pyro_rates(self, z=0, t=0, T0=float("inf"), betaKs=0.3333):
        """Compute the pyrolysis reaction rates at a given temperature"""
        temperature = T0 + betaKs * t
        k = self.generate_rate(temperature)
        dchidt = []
        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            if (1 - z[idx]) < 1e-5:
                z[idx] = 1
            dchidt.append(k[idx] * (1 - z[idx]) ** self.reaction_scheme_obj.dict_params["n"][idx])
        return dchidt

    def solve_system(self):
        """Solve the linear system of equation of pyrolysis reaction.
        By default, only Radau method for solving the initial value problem."""

        paramStep = 5  # for the max step in solver, adjust if needed
        max_step = paramStep / self.betaKs * 100
        temp_0 = self.temperature[0]

        y0 = np.zeros(self.reaction_scheme_obj.n_reactions)
        solution = solve_ivp(
            fun=lambda t, z: self.pyro_rates(z, t, temp_0, self.betaKs),
            t_span=(0, self.time[-1]),
            y0=y0,
            t_eval=self.time,
            method="Radau",
            max_step=max_step,
            rtol=1e-5,
        )

        # Compute the density evolution from chi
        percent_evo_sum = np.zeros(len(self.time))
        for chi, F in zip(solution.y, self.reaction_scheme_obj.dict_params["F"]):
            percent_evo = chi * F
            percent_evo_sum += percent_evo

        # Mass loss and mass loss rate
        self.rho_solid = self.reaction_scheme_obj.rhoIni * (1 - percent_evo_sum)
        self.drho_solid = np.gradient(-self.rho_solid, self.temperature)


@dataclass
class PyrolysisCompetitive(Pyrolysis):
    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False

    param_names: list = field(default_factory=list, init=False)
    dict_params: dict = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        # Convert beta in K/min to betaKs in K/sec
        self.betaKs = self.beta / 60
        # Compute temperature and time
        # Compute temperature and time
        self.time = compute_time(temp_0=self.temp_0, temp_end=self.temp_end, betaKs=self.betaKs, n_points=self.n_points)
        self.temperature = compute_temperature(self.time, temp_0=self.temp_0, betaKs=self.betaKs)

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def generate_matrix(self, temperature=float("inf")):
        """Build the coefficient matrix A for the linear pyrolysis reactions.
        See Eq. (3.6) in Torres et al., NASA TM)"""

        k_loss = np.zeros((self.reaction_scheme_obj.n_solids, self.reaction_scheme_obj.n_solids))
        k_gain = np.zeros((self.reaction_scheme_obj.n_solids, self.reaction_scheme_obj.n_solids))

        for solid in range(0, self.reaction_scheme_obj.n_solids):
            for reaction in range(0, self.reaction_scheme_obj.n_reactions):

                # Compute the loss in reactants
                if self.reaction_scheme_obj.solids[solid] == self.reaction_scheme_obj.solid_reactant[reaction]:
                    k_loss[solid][solid] += (10 ** self.reaction_scheme_obj.dict_params["A"][reaction]) * np.exp(
                        -self.reaction_scheme_obj.dict_params["E"][reaction] / (R * temperature)
                    )

                # Compute the gain in (solid) products
                if self.reaction_scheme_obj.solids[solid] == self.reaction_scheme_obj.solid_product[reaction]:
                    idx_reactant = self.reaction_scheme_obj.solids.index(
                        self.reaction_scheme_obj.solid_reactant[reaction]
                    )  # finds which reactant produces this solid
                    k_gain[solid][idx_reactant] += (
                            self.reaction_scheme_obj.g_sol[reaction]
                            * (10 ** self.reaction_scheme_obj.dict_params["A"][reaction])
                            * np.exp(-self.reaction_scheme_obj.dict_params["E"][reaction] / (R * temperature))
                    )  # k of Arrhenius
        return -k_loss + k_gain

    def pyro_rates(self, z, t):
        """Compute the pyrolysis reaction rates at a given temperature"""
        temperature = compute_temperature(t, temp_0=self.temp_0, betaKs=self.betaKs)
        drhodt = np.dot(self.generate_matrix(temperature), z)
        return drhodt

    def solve_system(self):
        """Solve the linear system of equation of pyrolysis reaction.
        By default, only Radau method for solving the initial value problem."""

        paramStep = 5  # for the max step in solver, adjust if needed
        max_step = paramStep / self.betaKs# * 100
        temp_0 = self.temperature[0]

        # Solve the system
        if self.isothermal:
            solution = solve_ivp(
                fun=lambda t, z: self.pyro_rates(z, t),
                t_span=(0, self.time[-1]),
                y0=self.reaction_scheme_obj.rhoIni,
                t_eval=self.time,
                rtol=1e-5,
            )
        else:
            solution = solve_ivp(
                fun=lambda t, z: self.pyro_rates(z, t),
                t_span=(0, self.time[-1]),
                y0=self.reaction_scheme_obj.rhoIni,
                t_eval=self.time,
                method="Radau",
                max_step=max_step,
                rtol=1e-5,
            )

        # Compute the total solid density
        for rho in solution.y:
            self.rho_solid += rho

        # # Check that dimensions are consistent
        # if len(self.rho_solid) != len(self.time):
        #     self.rho_solid = np.zeros(len(self.time))
        #     print("Inconsistency in rho_solid and self.time dimensions. rho_solid set to zero.")

        # Compute total solid density gradient
        self.drho_solid = np.gradient(-self.rho_solid, self.temperature)

        return solution.y
