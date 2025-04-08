from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import expi as ei

from pyropy.constants import R_gas  # exponential integral function

from .reaction_reader_writer import ReactManager


@dataclass
class Pyrolysis(Protocol):
    """
    Protocol for a Pyrolysis simulation class.

    Implementing classes should simulate pyrolysis processes and provide
    results such as solid density, temperature, and rate of change of solid density.
    """

    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False
    # for initialization
    rho_solid: np.ndarray = field(default_factory=lambda: np.empty(0))
    drho_solid: np.ndarray = field(default_factory=lambda: np.empty(0))
    time: np.ndarray = field(default_factory=lambda: np.empty(0))
    temperature: np.ndarray = field(default_factory=lambda: np.empty(0))

    def pyro_rates(self) -> list:
        """
        Compute the pyrolysis rates.

        Returns
        -------
        list of float
            A list of pyrolysis rates.
        """
        ...

    def solve_system(self) -> None:
        """
        Solve the system of equations for the pyrolysis process.

        This method should update `rho_solid`, `drho_solid`, `time`, and `temperature`.
        """
        ...

    def to_csv(self, filename: str) -> None:
        """
        Save the results to a CSV file.

        Parameters
        ----------
        filename : str
            The name of the CSV file to save the results to.

        Saves columns: 'time', 'temperature', 'rho', 'dRho'.
        """
        file_save = pd.DataFrame(
            {
                "time": self.time,
                "temperature": self.temperature,
                "rho": self.rho_solid,
                "dRho": self.drho_solid,
            }
        )
        file_save.to_csv(filename, index=False)


def compute_time(
    betaKs: float, temp_0: float = 373, temp_end: float = 2000, n_points: int = 500
) -> np.ndarray:
    """
    Compute a time vector for a linear temperature ramp.

    Parameters
    ----------
    betaKs : float
        Heating rate [K/s].
    temp_0 : float, optional
        Initial temperature [K]. Default is 373 K.
    temp_end : float, optional
        Final temperature [K]. Default is 2000 K.
    n_points : int, optional
        Number of time steps. Default is 500.

    Returns
    -------
    np.ndarray
        Time vector corresponding to the temperature ramp.
    """
    return np.linspace(0, (temp_end - temp_0) / betaKs, n_points)


def compute_temperature(time: np.ndarray, temp_0: float, betaKs: float) -> np.ndarray:
    """
    Compute temperature as a function of time for a linear heating ramp.

    Parameters
    ----------
    time : np.ndarray
        Time vector [s].
    temp_0 : float
        Initial temperature [K].
    betaKs : float
        Heating rate [K/s].

    Returns
    -------
    np.ndarray
        Temperature vector [K] over time.
    """
    return time * betaKs + temp_0


@dataclass
class PyrolysisParallel(Pyrolysis):
    """
    Pyrolysis model assuming parallel reactions.
    """

    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False

    param_names: list[str] = field(default_factory=list, init=False)
    dict_params: dict = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        # Convert beta in K/min to betaKs in K/sec
        self.betaKs = self.beta / 60
        # Compute temperature and time
        self.time = compute_time(
            temp_0=self.temp_0,
            temp_end=self.temp_end,
            betaKs=self.betaKs,
            n_points=self.n_points,
        )
        self.temperature = compute_temperature(
            self.time, temp_0=self.temp_0, betaKs=self.betaKs
        )

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def generate_rate(self, temperature: float = float("inf")) -> list[float]:
        """
        Compute reaction rates for all reactions at a given temperature.

        Parameters
        ----------
        temperature : float, optional
        Temperature [K] at which to evaluate the reaction rates.
        Default is infinity (i.e., returns zero rates due to exp(-inf)).

        Returns
        -------
        list[float]
        List of reaction rates.
        """
        k = []
        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            k.append(
                10 ** self.reaction_scheme_obj.dict_params["A"][idx]
                * np.exp(
                    -self.reaction_scheme_obj.dict_params["E"][idx]
                    / (R_gas * temperature)
                )
            )
        return k

    def pyro_rates(
        self,
        z: list[float] | None = None,
        t: float = 0.0,
        T0: float = float("inf"),
        betaKs: float = 0.3333,
    ) -> list[float]:
        """
        Compute pyrolysis reaction rates at a given time and base temperature.

        Parameters
        ----------
        z : list of float, optional
            List of conversion degrees for each reaction. If None, initializes to all zeros.
        t : float, optional
            Time [s] at which to compute the rates. Default is 0.0.
        T0 : float, optional
            Initial temperature [K]. Default is infinity.
        betaKs : float, optional
            Heating rate in [K/s]. Default is 0.3333.

        Returns
        -------
        List[float]
            Reaction rates (dα/dt) for each reaction.
        """
        if z is None:
            z = [0.0] * self.reaction_scheme_obj.n_reactions

        temperature = T0 + betaKs * t
        k = self.generate_rate(temperature)

        dchidt = []
        for idx in range(self.reaction_scheme_obj.n_reactions):
            z_val = min(z[idx], 1.0)  # Clip to avoid numerical errors
            exponent = self.reaction_scheme_obj.dict_params["n"][idx]
            dchidt.append(k[idx] * (1 - z_val) ** exponent)

        return dchidt

    def solve_system(self) -> None:
        """
        Solve the pyrolysis reaction system using the Radau method.

        Solves the system of ODEs for the conversion degrees of each reaction and
        computes the solid density and its time derivative.

        Notes
        -----
        The solver uses the Radau method, suitable for stiff systems.
        The temperature and time vectors must be initialized prior to calling this.
        """
        param_step = 5  # max step in K
        max_step = param_step / self.betaKs * 100
        temp_0 = self.temperature[0]
        y0 = np.zeros(self.reaction_scheme_obj.n_reactions)

        sol = solve_ivp(
            fun=lambda t, z: self.pyro_rates(z, t, temp_0, self.betaKs),
            t_span=(0, self.time[-1]),
            y0=y0,
            t_eval=self.time,
            method="Radau",
            max_step=max_step,
            rtol=1e-5,
        )

        percent_evo_sum = np.sum(
            [chi * F for chi, F in zip(sol.y, self.reaction_scheme_obj.dict_params["F"])],
            axis=0,
        )

        self.rho_solid = self.reaction_scheme_obj.rhoIni * (1 - percent_evo_sum)
        self.drho_solid = np.gradient(-self.rho_solid, self.temperature)


@dataclass
class PyrolysisCompetitive(Pyrolysis):
    """
    Pyrolysis model with competitive reactions.
    """

    temp_0: float
    temp_end: float
    beta: float
    n_points: int
    reaction_scheme_obj: ReactManager
    isothermal: bool = False

    param_names: list[str] = field(default_factory=list, init=False)
    dict_params: dict[str, npt.ArrayLike] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        # Convert beta in K/min to betaKs in K/sec
        self.betaKs = self.beta / 60
        # Compute temperature and time
        # Compute temperature and time
        self.time = compute_time(
            temp_0=self.temp_0,
            temp_end=self.temp_end,
            betaKs=self.betaKs,
            n_points=self.n_points,
        )
        self.temperature = compute_temperature(
            self.time, temp_0=self.temp_0, betaKs=self.betaKs
        )

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def generate_matrix(self, temperature: float = float("inf")) -> np.ndarray:
        """
        Build the coefficient matrix A for linear pyrolysis reactions.
        See Eq. (3.6) in Torres et al., NASA TM.

        Parameters
        ----------
        temperature : float
            Temperature at which the rate coefficients are evaluated (K).

        Returns
        -------
        np.ndarray
            Coefficient matrix A where A[i][j] describes the net rate of change of solid i
            due to reactions involving solid j.
        """

        n_solids = self.reaction_scheme_obj.n_solids
        n_reactions = self.reaction_scheme_obj.n_reactions

        k_loss = np.zeros((n_solids, n_solids))
        k_gain = np.zeros((n_solids, n_solids))

        for solid in range(n_solids):
            for reaction in range(n_reactions):

                A = 10 ** self.reaction_scheme_obj.dict_params["A"][reaction]
                E = self.reaction_scheme_obj.dict_params["E"][reaction]
                k = A * np.exp(-E / (R_gas * temperature))  # Arrhenius rate constant

                # Loss from this solid due to it being a reactant
                if (
                    self.reaction_scheme_obj.solids[solid]
                    == self.reaction_scheme_obj.solid_reactant[reaction]
                ):
                    k_loss[solid][solid] += k

                # Gain in this solid due to it being a product
                if (
                    self.reaction_scheme_obj.solids[solid]
                    == self.reaction_scheme_obj.solid_product[reaction]
                ):
                    idx_reactant = self.reaction_scheme_obj.solids.index(
                        self.reaction_scheme_obj.solid_reactant[reaction]
                    )
                    g = self.reaction_scheme_obj.g_sol[reaction]
                    k_gain[solid][idx_reactant] += g * k

        return -k_loss + k_gain

    def pyro_rates(self, z: np.ndarray, t: float) -> np.ndarray:
        """
        Compute the pyrolysis reaction rates using the system matrix A(T).

        Parameters
        ----------
        z : np.ndarray
            Vector of solid mass fractions or densities.
        t : float
            Time at which the rates are evaluated.

        Returns
        -------
        np.ndarray
            Time derivative of z, i.e., dz/dt at time t.
        """
        temperature = compute_temperature(t, temp_0=self.temp_0, betaKs=self.betaKs)  # type: ignore
        A_matrix = self.generate_matrix(temperature)  # type: ignore
        dzdt = np.dot(A_matrix, z)
        return dzdt

    def solve_system(self):
        """
        Solve the linear system of equations for pyrolysis reactions.
        Uses Radau method for non-isothermal conditions.
        Returns the full solution of solid densities over time.
        """

        paramStep = 5  # For max step size control
        max_step = paramStep / self.betaKs

        # Initial condition: vector of initial densities for each solid component
        y0 = np.array(self.reaction_scheme_obj.rhoIni)

        # Solve the system
        solution = solve_ivp(
            fun=lambda t, z: self.pyro_rates(z, t),
            t_span=(0, self.time[-1]),
            y0=y0,
            t_eval=self.time,
            method=(
                "Radau" if not self.isothermal else "RK45"
            ),  # RK45 is fine for isothermal if no stiffness
            max_step=max_step if not self.isothermal else np.inf,
            rtol=1e-5,
        )

        # Store each solid's individual solution if needed
        self.individual_rhos = (
            solution.y
        )  # Optional: if you want access to individual species

        # Total solid density (sum over all solids at each time point)
        self.rho_solid = np.sum(solution.y, axis=0)

        # Derivative of solid density with respect to temperature
        self.drho_solid = np.gradient(-self.rho_solid, self.temperature)

        return solution.y


@dataclass
class PyrolysisParallelAnalytical(Pyrolysis):
    """
    Analytical version of parallel reactions. See paper from Coheur.
    """

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
        self.time = compute_time(
            temp_0=self.temp_0,
            temp_end=self.temp_end,
            betaKs=self.betaKs,
            n_points=self.n_points,
        )
        self.temperature = compute_temperature(
            self.time, temp_0=self.temp_0, betaKs=self.betaKs
        )

        # initizalize rho and drho_solid
        self.rho_solid = np.zeros(self.n_points)  # initial density
        self.drho_solid = np.zeros(self.n_points)  # initial derivative density

    def pyro_rates(self): ...

    def solve_system(self):
        """For parallel reactions, there exists an analytical solution (see Torres, Coheur, NASA TM 2018).
        The solution is implemented here."""
        tau = self.betaKs
        T_0 = self.temperature[0]

        # Initialize arrays for cumulative evolution and reaction rate
        percent_evo_sum = np.zeros(len(self.time))
        pi_j = np.zeros(len(self.time))

        for idx in range(0, self.reaction_scheme_obj.n_reactions):
            xi_init = 0
            n = self.reaction_scheme_obj.dict_params["n"][idx]

            if n == 1:
                raise ZeroDivisionError("n cannot be equal to 1 in this mechanism!")

            A = 10 ** self.reaction_scheme_obj.dict_params["A"][idx]
            E = self.reaction_scheme_obj.dict_params["E"][idx]

            C = (
                (1 - xi_init) ** (1 - n) / (1 - n)
                + (A / tau) * T_0 * np.exp(-E / (R_gas * T_0))
                + ei(-E / (R_gas * T_0)) * E * (A / tau) / R_gas
            )

            xi_T = 1 - (
                (1 - n)
                * (
                    -(A / tau)
                    * self.temperature
                    * np.exp(-E / (R_gas * self.temperature))
                    - ei(-E / (R_gas * self.temperature)) * E * (A / tau) / R_gas
                    + C
                )
            ) ** (1 / (1 - n))

            # Update the cumulative evolution and reaction rate
            percent_evo_sum += xi_T * self.reaction_scheme_obj.dict_params["F"][idx]
            pi_j += (
                self.reaction_scheme_obj.dict_params["F"][idx]
                * (1 - xi_T) ** n
                * (A / tau)
                * np.exp(-E / (R_gas * self.temperature))
            )

        # Compute mass loss and mass loss rate
        self.rho_solid = self.reaction_scheme_obj.rhoIni * (
            1 - percent_evo_sum
        )  # Solid density evolution
        self.drho_solid = (
            self.reaction_scheme_obj.rhoIni * pi_j
        )  # Rate of change of solid density

        return self.rho_solid, self.drho_solid
