# import modules: SNowflake and OperatingConditions
from ethz_snow.snowfall import Snowfall
from ethz_snow.operatingConditions import OperatingConditions

# define the heat transfer parameters dictionary
d = {"int": 10, "ext": 100, "s0": 50, "s_sigma_rel": 0}

# define the cooling rate profile and holding steps
c = {"rate": 0.5 / 60, "start": 20, "end": -50}
h = []

# create an instance of the operating conditions class
op = OperatingConditions(t_tot=4 * 3600, cooling=c, holding=h)

# run a single simulation
S = Snowfall(
    k=d, opcond=op, seed_v=0
)  # seed_v is used to seed the rng for vial-dependent nuc parameters (See Snowflake docs)
S.run()

print(S.nucleationTimes())
print(S.nucleationTemperatures())
print(S.solidificationTimes())