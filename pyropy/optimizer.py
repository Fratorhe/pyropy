import os
import shutil
from typing import Type, Callable

import spotpy

from pyropy import ExperimentReader, ExperimentReaderCSV, Pyrolysis, PyrolysisParallel, ReactManager
from pyropy.auxiliary_functions import write_file_scheme, get_numbers_from_filename
from pyropy.rmse_multiple_files import rmse_multiple_files


class SpotpySetup(object):
    """Class to setup spotpy optimization routine

    Attributes:
        :atrr betas:
        :atrr params:
        :atrr names:
        :atrr times:
        :atrr dRho:
        :atrr Rho:
        :atrr temperatures:
        :atrr iternumber:
    """

    def __init__(
        self,
        files,
        params,
        folder,
        scheme_file,
        pyro_type: Type[Pyrolysis] = PyrolysisParallel,
        keepFolders=False,
        experiment_reader: Type[ExperimentReader] = ExperimentReaderCSV,
        isothermal=False,
        objective_function=Callable,
    ):
        """

        :param files: list of files to be treated
        :param params: list of params to be optimized
        :param folder: str folder where to read the experiments
        :param scheme_file: str with the scheme used
        :param pyro_type: str for type of pyrolysis (parallel, competitive, etc)
        :param keepFolders: bool if folder of simulations are saved
        :param isothermal: bool isothermal tests (not used in most cases)
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

        :param params:
        """
        for param in params:
            self.names.append(param.name)
        pass

    def parameters(self):
        """

        :return:
        """
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        """
        Simulate the case with ode_solver
        :param vector: vector of unknowns
        :return: list with [resultsdRho,resultsRho] from ode_solver
        """
        simulations = self.ode_solver(vector=vector)
        return simulations

    def evaluation(self):
        """

        :rtype: list with experimental observations
        """
        observations = [self.dRho, self.Rho]
        # observations = self.dRho
        return observations

    def objectivefunction(self, simulation, evaluation):
        """

        :param simulation: list with simulation results
        :param evaluation: list with experimental observations
        :return: float with value of the objective function
        """
        # Drho = rmse_multiple_files(evaluation,simulation) # 100 to give weight wrt rho
        rmse = rmse_multiple_files(evaluation, simulation)  # 100 to give weight wrt rho
        return rmse

    def ode_solver(self, vector):
        """

        :param vector: vector of unknowns
        :return: list of results for simulation function
        """
        results_dRho = []
        results_Rho = []
        self.iternumber += 1
        os.makedirs(str(self.iternumber))

        for beta, temperature in zip(self.betas, self.temperatures):
            write_file_scheme(
                filename=self.scheme_file, vector=vector, param_names=self.names, folder=str(self.iternumber) + "/"
            )
            reactions = ReactManager(filename=self.scheme_file, folder=str(self.iternumber) + "/")
            reactions.react_reader()
            reactions.param_reader()
            temperature = list(temperature)
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

        if self.keepFolders is False:
            shutil.rmtree(str(self.iternumber))

        return [results_dRho, results_Rho]
