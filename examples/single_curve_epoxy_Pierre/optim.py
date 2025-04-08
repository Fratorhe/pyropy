import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spotpy

from pyropy.auxiliary_functions import replace_results
from pyropy.experiment_reader import ExperimentReaderCSV
from pyropy.optimizer import SpotpySetup
from pyropy.pyrolysis import PyrolysisParallel
from pyropy.rmse_multiple_files import rmse_multiple_files

results = []
folder = "."
filesCSV = ["data/pyro_epoxy_rate_5.csv"]
params = [
    spotpy.parameter.Uniform("E1", low=1000, high=1e6, optguess=9.94796910e04),
    spotpy.parameter.Uniform("A1", low=1, high=15, optguess=4),
    spotpy.parameter.Uniform("n1", low=2, high=10, optguess=0.5),
    spotpy.parameter.Uniform("g1", low=0, high=1, optguess=0.3),
    spotpy.parameter.Uniform("Kcat1", low=0, high=1, optguess=0.3),
    spotpy.parameter.Uniform("F1", low=0, high=1, optguess=0.3),
]

spotpy_setup = SpotpySetup(
    files=filesCSV,
    params=params,
    folder=folder,
    scheme_file="reaction_scheme_optim.json",
    pyro_type=PyrolysisParallel,
    keepFolders=False,
    experiment_reader=ExperimentReaderCSV,
    objective_function=rmse_multiple_files,
)
rep = 500

sampler = spotpy.algorithms.sceua(
    spotpy_setup, dbname="SCEUA", dbformat="csv", save_sim=False, db_precision=np.float64
)
sampler.sample(rep)
results.append(sampler.getdata())
names = spotpy.analyser.get_parameternames(sampler.getdata())
best = spotpy.analyser.get_best_parameterset(sampler.getdata(), maximize=False)

# spotpy_setup.plotBest(bestResultVector=best[0])  #Doesn't work becuase 'spotpy_setup' object has no attribute 'plotBest'
bestResults = [x for x in best[0]]

resultCSV = pd.DataFrame.from_dict({"Params": names, "Value": bestResults})
resultCSV.to_csv(folder + "/" + "results_1")

replace_results(
    bestResults, names, "reaction_scheme_optim.json.template", "data_optimized.json"
)
