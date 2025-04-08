import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pyropy.pyrolysis import PyrolysisCompetitive
from pyropy.reaction_reader_writer import ReactManager


def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)


# stylePlots = style()
# plt.style.use(stylePlots)


figRho, axRho = plt.subplots()
figdRho, axdRho = plt.subplots()
figTime, axTime = plt.subplots()
figProduction, axProduction = plt.subplots()

axRho.set_xlabel("Temperature (K)")
axRho.set_ylabel("$\\rho/\\rho_{v} (-)$")
axdRho.set_xlabel("Temperature (K)")
axdRho.set_ylabel("d$(\\rho/\\rho_{v})$/d$T$ (mK$^{-1}$)")
axTime.set_xlabel("Temperature (K)")
axTime.set_ylabel("Time (s)")
axProduction.set_xlabel("Temperature (K)")
axProduction.set_ylabel("$\\rho/\\rho_{v} (-)$")

betas = (0.5, 5, 50, 500)

lineStyles = ["-", "--", "-.", ":"]

for rate, line in zip(betas, lineStyles):
    reactions = ReactManager("data_competing.json")

    reactions.react_reader()
    reactions.param_reader()
    test = PyrolysisCompetitive(
        temp_0=373,
        temp_end=1000,
        beta=rate,
        n_points=500,
        reaction_scheme_obj=reactions,
    )

    test.solve_system()

    rho = test.rho_solid
    drho = test.drho_solid
    t = test.time
    T = test.temperature

    axRho.plot(T, rho / 100, label=str(rate) + " K/min", linestyle=line, color="k")
    axdRho.plot(T, drho * 10, label=str(rate) + " K/min", linestyle=line, color="k")
    axTime.plot(T, t, label=str(rate) + " K/min", linestyle=line, color="k")

    # labels = ('Reactant', '$\\textrm{Activation}_{\\textrm{slow}}$', '$\\textrm{r}_{\\textrm{fast}}$','Solid 1','Solid 2')
    products = test.individual_rhos
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

figRho.savefig("rho.pdf")
figdRho.savefig("dRho.pdf")
figTime.savefig("Time.pdf")
figProduction.savefig("products.pdf")


plt.show()
# test.to_csv('test.csv')
