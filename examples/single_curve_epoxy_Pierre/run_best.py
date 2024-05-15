import pyropy.pyrolysis_epoxy as pyro
from pyropy.pyrolysis import PyrolysisParallel
import matplotlib.pyplot as plt

# from pybit import post_process

import pandas as pd

rates = [5]
# rates = [5,10,20,50,100]
T_0 = 308
T_end = 1373
n_T_steps = 100
for rate in rates:
    test = pyro.PyrolysisEpoxy(temp_0=T_0, temp_end=T_end, time=None, beta=rate, n_points=n_T_steps)

    test.react_reader(filename="data/data_epoxy.json")
    test.param_reader(filename="data/data_epoxy.json")
    test.solve_system()
    time = test.get_time()
    T = test.get_temperature()
    rho = test.get_rho_solid()
    drho = test.get_drho_solid()
    std_drho = max(drho) / 20

    plt.figure(1)
    plt.plot(T, rho)

    plt.figure(2)
    plt.plot(T, drho)
    # post_process.error_bar (T, drho, std_drho, 'C0')

    # Initialize pyrolysis model
    pyro_model = PyrolysisParallel(temp_0=T_0, temp_end=T_end, time=None, beta=rate, n_points=n_T_steps)

    # Read the parameters from the temporary file
    pyro_model.react_reader("data_optimized.json")
    pyro_model.param_reader("data_optimized.json")

    pyro_model.solve_system()
    time = pyro_model.get_time()
    T = pyro_model.get_temperature()
    rho = pyro_model.get_rho_solid()
    drho = pyro_model.get_drho_solid()
    std_drho = max(drho) / 20

    plt.figure(1)
    plt.plot(T, rho, "C1")

    plt.figure(2)
    plt.plot(T, drho, "C1")
    # post_process.error_bar (T, drho, std_drho, 'C1')

plt.show()

# Init guess {"E": 9.94796910e+04, "A": 9.42793325e+04,  "n": 0.5, "F":0.3,  "g":[1]}
