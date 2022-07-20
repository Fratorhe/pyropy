import matplotlib.pyplot as plt
import matplotlib.cm as colormaps

cmap = plt.cm.Dark2.colors

from pyropy import ReactManager, PyrolysisParallel, PyrolysisParallelAnalytical

reactions = ReactManager("data_parallel_verification.json")
reactions.react_reader()
reactions.param_reader()

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

beta = 5
test = PyrolysisParallel(temp_0=273, temp_end=2000, beta=beta, n_points=200, reaction_scheme_obj=reactions)

# Numerical solution using 200 points
test.solve_system()
rho = test.rho_solid
drho = test.drho_solid
t = test.time
T = test.temperature
test.to_csv(f"test_{beta}Kmin.csv")
ax1.plot(T, rho, label="Numerical")
ax2.plot(T, drho, label="Numerical")

# Analytical solution
test = PyrolysisParallelAnalytical(
    temp_0=273, temp_end=2000, beta=beta, n_points=200, reaction_scheme_obj=reactions
)
test.solve_system()
rho = test.rho_solid
drho = test.drho_solid
t = test.time
T = test.temperature
ax1.plot(T, rho, label="Analytical")
ax2.plot(T, drho, label="Analytical")


# Create the legend
ax1.set_title("Mass loss")
ax1.legend()

ax2.set_title("Mass loss rate")
ax2.legend()


plt.show()
