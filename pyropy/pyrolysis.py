from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
from scipy.integrate import solve_ivp

from .reaction_reader_writer import ReactManager

R = 8.314


class Pyrolysis(Protocol):
    def generate_rate(self) -> list:
        ...

    def pyro_rates(self) -> list:
        ...

    def solve_system(self) -> None:
        ...


def compute_time_temperature(betaKs: float, temp_0: float = 373, temp_end: float = 2000, n_points: int = 500):
    """ Compute temperature as a function of time based on the given heating rate and
    initial temperature temp_0. If no time vector is provided, a time vector of size n_points
    is computed from temp_0 to temp_end"""

    time = np.linspace(0, (temp_end - temp_0) / betaKs, n_points)
    temperature = time * betaKs + temp_0
    return time, temperature


@dataclass
class PyrolysisParallel:
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
        self.time, self.temperature = compute_time_temperature(temp_0=self.temp_0, temp_end=self.temp_end,
                                                               betaKs=self.betaKs, n_points=self.n_points)

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def generate_rate(self, temperature=float("inf")):
        k = []
        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            k.append(10 ** self.reaction_scheme_obj.dict_params['A'][idx] *
                     np.exp(-self.reaction_scheme_obj.dict_params['E'][idx] / (R * temperature)))
        return k

    def pyro_rates(self, z=0, t=0, T0=float("inf"), betaKs=0.3333):
        """ Compute the pyrolysis reaction rates at a given temperature """
        temperature = T0 + betaKs * t
        k = self.generate_rate(temperature)
        dchidt = []
        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            if (1 - z[idx]) < 1e-5:
                z[idx] = 1
            dchidt.append(k[idx] * (1 - z[idx]) ** self.reaction_scheme_obj.dict_params['n'][idx])
        return dchidt

    def solve_system(self):
        """ Solve the linear system of equation of pyrolysis reaction.
              By default, only Radau method for solving the initial value problem. """

        paramStep = 5  # for the max step in solver, adjust if needed
        max_step = paramStep / self.betaKs * 100
        temp_0 = self.temperature[0]

        y0 = np.zeros(self.reaction_scheme_obj.n_reactions)
        solution = solve_ivp(fun=lambda t, z: self.pyro_rates(z, t, temp_0, self.betaKs),
                           t_span=(0, self.time[-1]),
                           y0=y0,
                           t_eval=self.time,
                           method="Radau", max_step=max_step, rtol=1E-5)

        # Compute the density evolution from chi
        percent_evo_sum = np.zeros(len(self.time))
        for chi, F in zip(solution.y, self.reaction_scheme_obj.dict_params["F"]):
            percent_evo = chi * F
            percent_evo_sum += percent_evo

        # Mass loss and mass loss rate
        self.rho_solid = self.reaction_scheme_obj.rhoIni * (1 - percent_evo_sum)
        self.drho_solid = np.gradient(-self.rho_solid, self.temperature)
