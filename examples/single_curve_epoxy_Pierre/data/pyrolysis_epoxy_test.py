import pyrolysis_general.src.pyrolysis_epoxy as pyro
import matplotlib.pyplot as plt
from pybit import post_process

import pandas as pd

rates = [5]
# rates = [5,10,20,50,100]
for rate in rates:
    test = pyro.PyrolysisEpoxy(temp_0=308, temp_end=1373, time=None, beta=rate, n_points=100)

    test.react_reader(filename="data_epoxy.json")
    test.param_reader(filename="data_epoxy.json")
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
    post_process.error_bar(T, drho, std_drho, "C0")

    # Create the data list to be put in csv
    data = {"time": time, "temperature": T, "rho": rho, "dRho": drho, "std_drho": std_drho}

    df = pd.DataFrame(data, columns=data.keys())

    df.to_csv("pyro_epoxy_rate_" + str(rate) + ".csv")

# plt.plot(test.temperature, test.z.y[0])
# plt.plot(test.temperature, test.z.y[1])
# plt.plot(test.temperature, test.z.y[1]+test.z.y[0])
plt.show()
