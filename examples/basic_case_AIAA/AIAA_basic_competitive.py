from pyropy.pyrolysis import PyrolysisCompetitive
import matplotlib.pyplot as plt
from plot_python_vki import style
from matplotlib.lines import Line2D


def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)


stylePlots = style()
plt.style.use(stylePlots)


figRho, axRho = plt.subplots()
figdRho, axdRho = plt.subplots()
figTime, axTime = plt.subplots()
figProduction, axProduction = plt.subplots()

axRho.set_xlabel("Temperature (K)")
axRho.set_ylabel("$\\rho/\\rho_{\\textrm{v}} (-)$")
axdRho.set_xlabel("Temperature (K)")
axdRho.set_ylabel("d$(\\rho/\\rho_{\\textrm{v}})$/d$T$ (mK$^{-1}$)")
axTime.set_xlabel("Temperature (K)")
axTime.set_ylabel("Time (s)")
axProduction.set_xlabel("Temperature (K)")
axProduction.set_ylabel("$\\rho/\\rho_{\\textrm{v}} (-)$")

betas = (0.5, 5, 50, 500)

lineStyles = ["-", "--", "-.", ":"]

for rate, line in zip(betas, lineStyles):
    test = PyrolysisCompetitive(temp_0=373, temp_end=1000, time=None, beta=rate, n_points=500)
    test.react_reader("data_competing.json")
    test.param_reader("data_competing.json")
    test.solve_system()

    rho = test.get_rho_solid()
    drho = test.get_drho_solid()
    t = test.get_time()
    T = test.get_temperature()

    axRho.plot(T, rho / 100, label=str(rate) + " K/min", linestyle=line, color="k")
    axdRho.plot(T, drho * 10, label=str(rate) + " K/min", linestyle=line, color="k")
    axTime.plot(T, t, label=str(rate) + " K/min", linestyle=line, color="k")

    # labels = ('Reactant', '$\\textrm{Activation}_{\\textrm{slow}}$', '$\\textrm{r}_{\\textrm{fast}}$','Solid 1','Solid 2')
    products = test.z.y
    axProduction.plot(T, products[1] / 100, color="C0", linestyle=line)
    axProduction.plot(T, products[2] / 100, color="C3", linestyle=line)

    # figProduction.savefig(str(rate)+'.eps')
axRho.legend(loc="best")
axdRho.legend(loc="best")
axTime.legend(loc="best")
axProduction.legend(loc="best")

# Create the legend
lines = [
    ("Solid$_1$", {"color": "C0", "linestyle": "-"}),
    ("Solid$_2$", {"color": "C3", "linestyle": "-"}),
    ("", {"color": "C3", "linestyle": "None"}),
    ("0.5 K/min", {"color": "k", "linestyle": "-"}),
    ("5 K/min", {"color": "k", "linestyle": "--"}),
    ("50 K/min", {"color": "k", "linestyle": "-."}),
    ("500 K/min", {"color": "k", "linestyle": ":"}),
]

legend1 = plt.legend(
    # Line handles
    [create_dummy_line(**l[1]) for l in lines],
    # Line titles
    [l[0] for l in lines],
    loc="best",
    frameon=False,
)

axProduction.add_artist(legend1)

for save_format in (".png", ".eps"):
    figRho.savefig("rho" + save_format)
    figdRho.savefig("dRho" + save_format)
    figTime.savefig("Time" + save_format)
    figProduction.savefig("products" + save_format)


plt.show()
# test.to_csv('test.csv')
