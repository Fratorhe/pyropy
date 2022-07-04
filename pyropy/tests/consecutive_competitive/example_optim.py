import numpy as np
import pandas as pd
import spotpy

from pyropy import SpotpySetup, PyrolysisCompetitive, ExperimentReaderCSV, rmse_multiple_files, replace_results

results = []
folder = "."

# filesCSV = ["test_20Kmin.csv", "test_1Kmin.csv"]
filesCSV = ["test_5Kmin.csv"]

params = [
    spotpy.parameter.Uniform('A1', low=1, high=7, optguess=3.9),
    spotpy.parameter.Uniform('E1', low=50E3, high=100E3, optguess=60E3),
    spotpy.parameter.Uniform('A2', low=1, high=7, optguess=3.9),
    spotpy.parameter.Uniform('E2', low=80E3, high=120E3, optguess=90E3),
]

spotpy_setup = SpotpySetup(files=filesCSV, params=params,
                           folder=folder, scheme_file="data_competitive_verification.json",
                           pyro_type=PyrolysisCompetitive, keepFolders=False, experiment_reader=ExperimentReaderCSV,
                           objective_function=rmse_multiple_files)

rep=600

sampler = spotpy.algorithms.sceua(spotpy_setup, dbname='SCEUA', dbformat='csv', save_sim=False, db_precision=np.float64)
sampler.sample(rep)
results.append(sampler.getdata())
names = spotpy.analyser.get_parameternames(sampler.getdata())
best = spotpy.analyser.get_best_parameterset(sampler.getdata(), maximize=False)

bestResults = [x for x in best[0]]

resultCSV = pd.DataFrame.from_dict({"Params":names, "Value":bestResults})
print(resultCSV)
resultCSV.to_csv(folder+'/'+"results_1")

replace_results(bestResults, names, 'data_competitive_verification.json.template', 'data_optimized.json')
