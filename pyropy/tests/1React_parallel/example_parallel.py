import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
cmap = plt.cm.Dark2.colors

from pyropy import ReactManager, PyrolysisParallel

reactions = ReactManager("data_parallel_verification.json")
reactions.react_reader()
reactions.param_reader()


beta = 5
test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=beta, n_points=200, reaction_scheme_obj=reactions)

# Numerical solution using 200 points
test.solve_system()
rho = test.rho_solid
drho = test.drho_solid
t = test.time
T = test.temperature
test.to_csv(f'test_{beta}Kmin.csv')
plt.figure(1)
plt.plot(T, rho,label='Numerical')
plt.figure(2)
plt.plot(T, drho,label='Numerical')

# Analytical solution
test.compute_analytical_solution()
rho = test.rho_solid
drho = test.drho_solid
t = test.time
T = test.temperature
plt.figure(1)
plt.plot(T, rho,label='Analytical')
plt.figure(2)
plt.plot(T, drho,label='Analytical')


# Create the legend
plt.figure(1)
plt.title('Mass loss')
plt.legend()

plt.figure(2)
plt.title('Mass loss rate')
plt.legend()


plt.show()
