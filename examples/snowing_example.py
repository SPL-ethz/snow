# import modules
from ethz_snow.snowing import Snowing
from ethz_snow.operatingConditions import OperatingConditions

# define the heat transfer parameters dictionary
d = {"int": 20, "ext": 0, "s0": 50, "s_sigma_rel": 0}

# define the cooling rate profile and holding steps
c = {"rate": 0.5 / 60, "start": 20, "end": -50}
h = []

# create an instance of the operating conditions class
op = OperatingConditions(t_tot=3 * 3600, cooling=c, holding=h)

# run a single spatial simulation
S = Snowing(k=d, opcond=op)
S.run()
