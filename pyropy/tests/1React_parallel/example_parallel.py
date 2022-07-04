import matplotlib.pyplot as plt

from pyropy import ReactManager, PyrolysisParallel

reactions = ReactManager("data_parallel_verification.json")
reactions.react_reader()
reactions.param_reader()


beta = 5
test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=beta, n_points=200, reaction_scheme_obj=reactions)


test.solve_system()

rho = test.rho_solid
t = test.time
T = test.temperature
plt.plot(T, rho)
plt.show()
test.to_csv(f'test_{beta}Kmin.csv')
