import matplotlib.pyplot as plt

from pyropy.pyrolysis import PyrolysisCompetitive
from pyropy.reaction_reader_writer import ReactManager

reactions = ReactManager("data_competing.json")
reactions.react_reader()
reactions.param_reader()

test = PyrolysisCompetitive(
    temp_0=373, temp_end=2000, beta=1, n_points=500, reaction_scheme_obj=reactions
)

# test.react_writer("tests")
test.solve_system()
# test.plot_solid_density()
#
rho = test.rho_solid
t = test.time
plt.plot(t, rho)
plt.show()
# test.to_csv('test.csv')
