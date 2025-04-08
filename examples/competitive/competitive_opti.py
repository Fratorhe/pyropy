import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spotpy

from pyropy.auxiliary_functions import replace_results
from pyropy.experiment_reader import ExperimentReaderCSV
from pyropy.optimizer import SpotpySetup
from pyropy.pyrolysis import PyrolysisCompetitive
from pyropy.rmse_multiple_files import rmse_multiple_files

results = []
folder = "."
# filesCSV = [f for f in os.listdir(folder) if f.endswith('.csv')]
filesCSV = ["test_20Kmin.csv"]
# filesCSV = ["Wong_10Kmin_new.csv"]
# filesCSV = ["Wong_10Kmin_new.csv", "Wong_10Kmin_new.csv"]
params = [
    spotpy.parameter.Uniform("E1", low=25e3, high=60e3, optguess=47000),
    spotpy.parameter.Uniform("E2", low=100e3, high=170e3, optguess=150000),
    spotpy.parameter.Uniform("E3", low=40e3, high=80e3, optguess=50000),
    spotpy.parameter.Uniform("E4", low=40e3, high=60e3, optguess=51000),
    spotpy.parameter.Uniform("A1", low=1, high=6, optguess=3.9),
    spotpy.parameter.Uniform("gamma1", low=0.1, high=0.2, optguess=0.99),
    spotpy.parameter.Uniform("gamma2", low=0.005, high=0.02, optguess=0.99),
]

spotpy_setup = SpotpySetup(
    files=filesCSV,
    params=params,
    folder=folder,
    scheme_file="data_competing.json",
    pyro_type=PyrolysisCompetitive,
    keepFolders=False,
    experiment_reader=ExperimentReaderCSV,
    objective_function=rmse_multiple_files,
)
rep = 200

sampler = spotpy.algorithms.sceua(
    spotpy_setup, dbname="SCEUA", dbformat="csv", save_sim=False, db_precision=np.float64
)
sampler.sample(rep)
results.append(sampler.getdata())
names = spotpy.analyser.get_parameternames(sampler.getdata())
best = spotpy.analyser.get_best_parameterset(sampler.getdata(), maximize=False)

bestResults = [x for x in best[0]]

resultCSV = pd.DataFrame.from_dict({"Params": names, "Value": bestResults})
resultCSV.to_csv(folder + "/" + "results_1")

replace_results(
    bestResults,
    names,
    "data_competing.json.template",
    "data_optimized.json",
)
