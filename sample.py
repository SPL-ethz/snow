# import all modules
from ethz_snow.snowflake import Snowflake
from ethz_snow.snowfall import Snowfall
from ethz_snow.snowing import Snowing
from ethz_snow.operatingConditions import OperatingConditions

# parameters and operating conditions
d = {"int": 20, "ext": 0, "s0": 20, "s_sigma_rel": 0.1}
c = {"rate": 0.5 / 60, "start": 20, "end": -50}
h = [dict(duration=60 * 60, temp=-10), dict(duration=60 * 60, temp=-5)]
op = OperatingConditions(t_tot=3e4, cooling=c, holding=h, cnTemp=-5)

# snowflake
S_flake = Snowflake(k=d, Nrep=10, N_vials=(7, 7, 1), opcond=op)
S_flake.run()

# snowfall, shelf-scale
S_shelf = Snowfall(pool_size=8, k=d, Nrep=10, N_vials=(7, 7, 1), opcond=op)
S_shelf.run()

# snowfall, pallet-scale
d = {"int": 10, "ext": 10, "s0": 0, "s_sigma_rel": 0}
c = {
    "rate": 1e-16,
    "start": -8,
    "end": -8,
}  # rate equals zero yields divide by zero error
initial = {"temp": 20}
op = OperatingConditions(t_tot=6e6, cooling=c, holding=dict(duration=6e6, temp=-8))
S_pallet = Snowfall(
    pool_size=8,
    k=d,
    Nrep=128,
    N_vials=(40, 36, 18),
    opcond=op,
    dt=5,
    initialStates=initial,
)
S_pallet.run()

# snowing
d = {"int": 0, "ext": 0, "s0": 20, "s_sigma_rel": 0}
c = {"rate": 1 / 60, "start": 20, "end": -50}
h = []
op = OperatingConditions(t_tot=3 * 3600, cooling=c, holding=h)
S_ing = Snowing(
    k=d, opcond=op, temperature="spatial_1D", configuration="shelf", plotting=True
)
S_ing.run()
