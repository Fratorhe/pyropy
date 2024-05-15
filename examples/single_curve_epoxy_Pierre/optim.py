import matplotlib.pyplot as plt
import spotpy
from pyropy.optimizer import spotpy_setup
import pandas as pd
import pyropy.pyrolysis as pyro

results = []
folder = "."
filesCSV = ["data/pyro_epoxy_rate_5.csv"]
params = [
    spotpy.parameter.Uniform("E1", low=1000, high=1e6, optguess=9.94796910e04),
    spotpy.parameter.Uniform("A1", low=1, high=15, optguess=4),
    spotpy.parameter.Uniform("n1", low=2, high=10, optguess=0.5),
    spotpy.parameter.Uniform("g1", low=0, high=1, optguess=0.3),
]

spotpy_setup = spotpy_setup(
    files=filesCSV,
    params=params,
    folder=folder,
    scheme_file="reaction_scheme_optim.json",
    pyro_type="PyrolysisParallel",
    keepFolders=False,
)
rep = 500

sampler = spotpy.algorithms.sceua(
    spotpy_setup, dbname="SCEUA", dbformat="ram", save_sim=False, alt_objfun="rmse_multiple_files"
)
sampler.sample(rep, ngs=4)
results.append(sampler.getdata())
names = spotpy.analyser.get_parameternames(sampler.getdata())
best = spotpy.analyser.get_best_parameterset(sampler.getdata())

# spotpy_setup.plotBest(bestResultVector=best[0])  #Doesn't work becuase 'spotpy_setup' object has no attribute 'plotBest'
bestResults = [x for x in best[0]]

resultCSV = pd.DataFrame.from_dict({"Params": names, "Value": bestResults})
resultCSV.to_csv(folder + "/" + "results_1")

pyro.replace_results(bestResults, names, "reaction_scheme_optim.json.template", "data_optimized.json")
