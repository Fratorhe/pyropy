import os
import shutil
from typing import Callable, Type

import spotpy

from pyropy import (
    ExperimentReader,
    ExperimentReaderCSV,
    PyrolysisParallel,
    ReactManager,
)
from pyropy.auxiliary_functions import get_numbers_from_filename, write_file_scheme
from pyropy.rmse_multiple_files import rmse_multiple_files


class SpotpySetup:
    """
    Class to set up and run SPOTPY optimization for pyrolysis models.

    Attributes
    ----------
    betas : list of float
        Heating rates extracted from filenames.
    params : list
        SPOTPY parameter definitions.
    names : list of str
        Names of parameters to optimize.
    times : list of np.ndarray
        Time arrays from experiments.
    dRho : list of np.ndarray
        Experimental dRho data.
    Rho : list of np.ndarray
        Experimental Rho data.
    temperatures : list of np.ndarray
        Experimental temperature profiles.
    iternumber : int
        Counter for simulation iterations.
    """

    def __init__(
            self,
            files: list[str],
            params: list,
            folder: str,
            scheme_file: str,
            pyro_type: Type = PyrolysisParallel,
            keepFolders: bool = False,
            experiment_reader: Type[ExperimentReader] = ExperimentReaderCSV,
            isothermal: bool = False,
            objective_function: Callable | None = None,
    ) -> None:
        """
        Initialize the SPOTPY setup.

        Parameters
        ----------
        files : list of str
            Filenames with experimental data.
        params : list
            Parameters to be optimized.
        folder : str
            Directory with experimental files.
        scheme_file : str
            Path to reaction scheme template.
        pyro_type : Type
            Pyrolysis model class (default: PyrolysisParallel).
        keepFolders : bool
            Whether to retain folders created during simulations.
        isothermal : bool
            Whether simulations are isothermal (default: False).
        objective_function : Callable
            Objective function to compare simulations with experiments.
        """
        self.pyro_type = pyro_type
        self.folder = folder
        self.files = files
        self.betas = []
        for filename in self.files:
            self.betas.append(float(get_numbers_from_filename(filename)[0]))
        self.params = params
        self.names = []
        self.get_param_names(params)
        self.times = []
        self.dRho = []
        self.Rho = []
        self.temperatures = []
        self.scheme_file = scheme_file

        for filename in self.files:
            experiment = experiment_reader(filename=filename, folder=folder)
            self.times.append(experiment.time)
            self.temperatures.append(experiment.temperature)
            self.dRho.append(experiment.dRho)
            self.Rho.append(experiment.Rho)

        self.iternumber = 0
        self.keepFolders = keepFolders
        self.pyro_type = pyro_type
        self.objective_function = objective_function
        # X = self.dRho+self.Rho

    def get_param_names(self, params):
        """
        Extract and store the names of parameters.

        Parameters
        ----------
        params : list
            List of parameter objects, where each object has a 'name' attribute.
        """
        for param in params:
            self.names.append(param.name)
        pass

    def parameters(self):
        """
        Generate a set of SPOTPY parameters.

        Returns
        -------
        spotpy.parameter.ParameterSet
            A set of parameters ready for SPOTPY optimization.
        """
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector: list[float]) -> list:
        """
        Run the simulation with the provided parameter vector.

        Parameters
        ----------
        vector : list of float
            Vector of parameter values to be used in the simulation.

        Returns
        -------
        list
            List containing simulation results [results_dRho, results_Rho].
        """
        simulations = self.ode_solver(vector=vector)
        return simulations

    def evaluation(self) -> list:
        """
        Retrieve experimental observations for comparison.

        Returns
        -------
        list
            A list containing the experimental data [dRho, Rho].
        """
        return [self.dRho, self.Rho]

    def objectivefunction(self, simulation: list, evaluation: list) -> float:
        """
        Calculate the objective function value (e.g., RMSE) between simulation and experimental data.

        Parameters
        ----------
        simulation : list
            List containing simulation results [dRho, Rho].
        evaluation : list
            List containing experimental data [dRho, Rho].

        Returns
        -------
        float
            The computed objective function value (e.g., RMSE).
        """
        return rmse_multiple_files(evaluation, simulation)

    def ode_solver(self, vector: list[float]) -> list[list]:
        """
        Solve the ODE system for each experiment using the provided parameter vector.

        Parameters
        ----------
        vector : list of float
            A list of parameter values (unknowns) for the simulation.

        Returns
        -------
        list
            A list containing two sublists: [results_dRho, results_Rho] from simulations.
        """
        results_dRho, results_Rho = [], []
        self.iternumber += 1
        os.makedirs(str(self.iternumber))

        for beta, temperature in zip(self.betas, self.temperatures):
            sim_folder = f"{self.iternumber}/"
            write_file_scheme(
                filename=self.scheme_file,
                vector=vector,
                param_names=self.names,
                folder=sim_folder,
            )
            reactions = ReactManager(filename=self.scheme_file, folder=sim_folder)
            reactions.react_reader()
            reactions.param_reader()

            n_timesteps = len(temperature)
            simulation = self.pyro_type(
                temp_0=temperature[0],
                temp_end=temperature[-1],
                beta=beta,
                n_points=n_timesteps,
                isothermal=False,
                reaction_scheme_obj=reactions,
            )
            simulation.solve_system()
            results_dRho.append(simulation.drho_solid)
            results_Rho.append(simulation.rho_solid)

        if not self.keepFolders:
            shutil.rmtree(str(self.iternumber))

        return [results_dRho, results_Rho]
