import matplotlib.pyplot as plt

from pyropy import ReactManager, PyrolysisCompetitive

reactions = ReactManager("data_competitive_verification.json")
reactions.react_reader()
reactions.param_reader()


beta = 5
test = PyrolysisCompetitive(temp_0=273, temp_end=1500, beta=beta, n_points=200, reaction_scheme_obj=reactions)


y = test.solve_system()

rho = test.rho_solid
t = test.time
T = test.temperature

plt.plot(T, rho, label="total")
plt.plot(T, y[0, :], label="0")
plt.plot(T, y[1, :], label="1")
plt.plot(T, y[2, :], label="2")
plt.legend()
plt.show()
test.to_csv(f"test_{beta}Kmin.csv")
