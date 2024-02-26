# import modules: Snowing and OperatingConditions
from ethz_snow.snowing import Snowing
from ethz_snow.operatingConditions import OperatingConditions

# define the heat transfer parameters dictionary
d = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}

# define the cooling rate profile and holding steps
c = {"rate": 0.5 / 60, "start": 20, "end": -50}
h = []

# create an instance of the operating conditions class
op = OperatingConditions(t_tot=4 * 3600, cooling=c, holding=h)


# run a single spatial simulation
S_1 = Snowing(k=d, opcond=op, Nrep=1)
S_1.run()

# show results
S_1.results

# plot evolutions of temperature and ice mass fraction
S_1.plot_evolution(what="temperature")
S_1.plot_evolution(what="ice_mass_fraction")

# get individual arrays
time = S_1.time
shelf = S_1.shelfTemp
temp = S_1.temp
ice = S_1.iceMassFraction


# multiple simulations are run if Nrep > 1
S = Snowing(k=d, opcond=op, Nrep=50)
S.run()
S.results

# plot
S.plot_cdf(what="T_nuc")
S.plot_cdf(what="t_nuc")
S.plot_cdf(what="t_sol")
S.plot_cdf(what="t_fr")
