import matplotlib.pyplot as plt

from pyropy.pyrolysis_epoxy import PyrolysisEpoxy
from pyropy.reaction_reader_writer import ReactManager

rates = [5]
# rates = [5,10,20,50,100]
T_0 = 308
T_end = 1373
n_T_steps = 100
reactions = ReactManager("data_optimized.json")
reactions.react_reader()
reactions.param_reader()

for rate in rates:
    test = PyrolysisEpoxy(
        temp_0=T_0,
        temp_end=T_end,
        beta=rate,
        n_points=n_T_steps,
        reaction_scheme_obj=reactions,
    )

    test.solve_system()
    time = test.time
    T = test.temperature
    rho = test.rho_solid
    drho = test.drho_solid
    std_drho = max(drho) / 20

    plt.figure(1)
    plt.plot(T, rho)

    plt.figure(2)
    plt.plot(T, drho)
    # post_process.error_bar (T, drho, std_drho, 'C0')


plt.show()

# Init guess {"E": 9.94796910e+04, "A": 9.42793325e+04,  "n": 0.5, "F":0.3,  "g":[1]}
