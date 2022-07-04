import matplotlib.pyplot as plt

from pyropy import PyrolysisParallel
from pyropy import ReactManager

reactions = ReactManager(filename="data_parallel_verification.json")
reactions.react_reader()
reactions.param_reader()
test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=1, n_points=500, reaction_scheme_obj=reactions)
test.solve_system()
#
rho = test.rho_solid
T = test.temperature
plt.plot(T, rho)

test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=5, n_points=500, reaction_scheme_obj=reactions)
test.solve_system()
#
rho = test.rho_solid
T = test.temperature
plt.plot(T, rho)

reactions = ReactManager(filename="data_optimized.json")
reactions.react_reader()
reactions.param_reader()
test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=1, n_points=500, reaction_scheme_obj=reactions)
test.solve_system()
rho = test.rho_solid
T = test.temperature
plt.plot(T, rho)

test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=5, n_points=500, reaction_scheme_obj=reactions)
test.solve_system()
rho = test.rho_solid
T = test.temperature
plt.plot(T, rho)

plt.show()
