# import modules: Snowflake and OperatingConditions
from ethz_snow.snowflake import Snowflake
from ethz_snow.operatingConditions import OperatingConditions

# define the heat transfer parameters dictionary
d = {"int": 20, "ext": 0, "s0": 50, "s_sigma_rel": 0}

# define the cooling rate profile and holding steps
c = {"rate": 0.5 / 60, "start": 20, "end": -50}
h = []

# create an instance of the operating conditions class
op = OperatingConditions(t_tot=3 * 3600, cooling=c, holding=h)

# run simulation
S = Snowflake(k=d, opcond=op, N_vials=(5, 5, 1))
S.run()