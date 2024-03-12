#  ____     _   _     ___    __        __   _
# / ___|   | \ | |   / _ \   \ \      / /  (_)   _ __      __ _
# \___ \   |  \| |  | | | |   \ \ /\ / /   | |  | '_ \    / _` |
#  ___) |  | |\  |  | |_| |    \ V  V /    | |  | | | |  | (_| |
# |____/   |_| \_|   \___/      \_/\_/     |_|  |_| |_|   \__, |
#                                                         |___/

"""Implementation of the SNOWing class.

This module contains the SNOWing class used to run the single vial simulation
with spatial information.
"""

import numpy as np
import pandas as pd
import seaborn as sns
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.integrate import simps
from typing import Optional
from scipy.stats import norm
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerTuple

import ethz_snow.utils as Utils
from ethz_snow.constants import calculateDerived
from ethz_snow.operatingConditions import OperatingConditions

HEATFLOW_REQUIREDKEYS = ("int", "ext", "s0")
VIAL_GROUPS = ("corner", "edge", "center", "all")


class Snowing:
    """A class to handle a single SNOWing (Stochastic Nucleation of Water
     with INternal Gradients) simulation.

    More information regarding the equations and their derivation can be found in
    https://doi.org/10.1016/j.cej.2024.148660, Deck et al. (2024).

    Parameters:
        k (dict): Heat transfer coefficients.
        opcond (OperatingConditions): Operating conditions of the run.
        temperature (str): Model dimensionality.
        configuration (str): Freezing configuration.
        plotting (bool): Whether evolutions are plotted or not.
        Nrep (int): Number of repetitions.
        configPath (Optional[str]): The path of the (optional) custom config yaml.
    """

    def __init__(
        self,
        k: dict = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0},
        opcond: OperatingConditions = OperatingConditions(),
        Nrep: int = 1,
        configPath: Optional[str] = None,
    ):
        """Construct a Snowing object.

        Args:
            k (dict, optional): Heat transfer coefficients.
                Defaults to {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}.
            opcond (OperatingConditions, optional):
                Operating conditions of the run. Defaults to OperatingConditions().
            Nrep (int, optional): Number of repetitions. Defaults to 1.
            configPath (Optional[str], optional): The path of the (optional)
                custom config yaml. Defaults to None.

        Raises:
            TypeError: If k is not a dict.
            ValueError: If a aspecific key is missing in dict k.
            TypeError: If opcond is not of type operatingConditions.

        Examples:
            The following command will create and instance of the Snowing class with
            default operating conditions and heat transfer parameters:

            >>> S = Snowing()

            To customize the operating conditions, first define the operating conditions
            and the heat transfer parameters:

            >>> d = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}
            >>> c = {"rate": 0.5 / 60, "start": 20, "end": -50}
            >>> h = [dict(duration=60*60, temp=-10)]
            >>> op = OperatingConditions(t_tot=5*3600, cooling=c, holding=h)

            Then an instance of the Snowing class can be created:

            >>> S = Snowing(k=d, opcond=op)

            See tutorial for more information on how to use the spatial model functionalities.
        """

        if not isinstance(k, dict):
            raise TypeError(f"Input k must be of type dict. Was given {type(k)}.")
        elif not all([key in k.keys() for key in HEATFLOW_REQUIREDKEYS]):
            raise ValueError(
                (
                    f"A required key was missing from dictionary k, specifically "
                    + f"{set(HEATFLOW_REQUIREDKEYS) - set(k.keys())}."
                )
            )
        else:
            self.k = k

        # for saving the results
        self._domain = None
        self._prop = None
        self._time = None
        self._shelfTemp = None
        self._temp = None
        self._iceMassFraction = None

        # dictionary of results
        self._stats = {}

        # number of repetitions
        self.Nrep = Nrep

        # dictionary of results for multiple simulations
        if self.Nrep > 1:
            self._statsMultiple = {}

        # path to the parameters file
        self.configPath = configPath

        # simulation progress
        self._simulationStatus = 0

        # operating conditions
        if not isinstance(opcond, OperatingConditions):
            raise TypeError(
                "Input opcond must be an instance of class OperatingConditions."
            )
        else:
            self.opcond = opcond

    @property
    def configPath(self) -> Optional[str]:
        """Get configPath property."""
        return self._configPath

    @configPath.setter
    def configPath(self, value):
        """Set configPath property."""
        self.const = calculateDerived(value)
        self._configPath = value

    @property
    def simulationStatus(self) -> int:
        """Return simulation status of instance.

        Returns:
            int: Simulation status. 0 = not run, 1 = run.
        """
        if (self._simulationStatus == 1) or (self._temp is not None):
            self._simulationStatus = 1

        return self._simulationStatus

    @property
    def time(self) -> np.ndarray:
        """Return time array in the simulation of instance.

        Returns:
            np.ndarray: Time array in hours.
        """
        assert (self._simulationStatus == 1) or (
            self._time is not None
        ), "Simulation not yet carried out or completed."
        return self._time

    @property
    def temp(self) -> np.ndarray:
        """Return temperature evolution in the simulation of instance.

        Returns:
            np.ndarray: Temperature evolution in °C.
        """
        assert (self._simulationStatus == 1) or (
            self._temp is not None
        ), "Simulation not yet carried out or completed."
        return self._temp

    @property
    def shelfTemp(self) -> np.ndarray:
        """Return shelf temperature evolution in the simulation of instance.

        Returns:
            np.ndarray: Shelf temperature evolution in °C.
        """
        assert (self._simulationStatus == 1) or (
            self._shelfTemp is not None
        ), "Simulation not yet carried out or completed."
        return self._shelfTemp

    @property
    def iceMassFraction(self) -> np.ndarray:
        """Return ice mass fraction evolution in simulation of instance.

        Returns:
            np.ndarray: Ice mass fraction evolution (no units).
        """
        assert (self._simulationStatus == 1) or (
            self._iceMassFraction is not None
        ), "Simulation not yet carried out or completed."
        return self._iceMassFraction

    @property
    def results(self) -> dict:
        """Return minimum, mean, kinetic mean and maximum nucleation temperature
            in simulation of instance.

        Returns:
            dict: 4 values of nucleation temperature (minimum, kinetic mean,
            mean and maximum). Irrelevant for homogeneous model simulation.
        """
        assert (
            self._simulationStatus == 1
        ), "Simulation not yet carried out or completed."
        if self.Nrep == 1:
            self._stats_df = pd.DataFrame.from_dict(self._stats, orient="index").T
            return self._stats_df
        elif self.Nrep > 1:
            self._statsMultiple_df = pd.DataFrame.from_dict(
                self._statsMultiple, orient="index"
            )
            self._statsMultiple_df.columns = list(self._keys)
            return self._statsMultiple_df

    # simulate the entire process
    def run(self, how="async"):
        # homogeneous model
        if self.const["dimensionality"] == "homogeneous":
            if self.Nrep == 1:
                self._run_0D()
            else:
                if how == "sequential":
                    for i in range(self.Nrep):
                        self._run_0D(seed=i)
                elif how == "async":
                    with mp.Pool(mp.cpu_count()) as pool:
                        res = pool.starmap_async(
                            self._run_0D,
                            [(i,) for i in range(self.Nrep)],
                        ).get()
                        self._statsMultiple = {i: r for i, r in enumerate(res)}
                self._keys = ["T_nuc", "t_nuc", "t_sol", "t_fr"]

        # 1D spatial model
        elif self.const["dimensionality"] == "spatial_1D":
            if self.Nrep == 1:
                self._run_1D()
            else:
                if how == "sequential":
                    for i in range(self.Nrep):
                        self._run_1D(seed=i)
                elif how == "async":
                    with mp.Pool(mp.cpu_count()) as pool:
                        res = pool.starmap_async(
                            self._run_1D,
                            [(i,) for i in range(self.Nrep)],
                        ).get()
                        self._statsMultiple = {i: r for i, r in enumerate(res)}
                self._keys = [
                    "T_nuc_min",
                    "T_nuc_kin",
                    "T_nuc_mean",
                    "T_nuc_max",
                    "t_nuc",
                    "t_sol",
                    "t_fr",
                ]

        # 2D spatial model
        elif self.const["dimensionality"] == "spatial_2D":
            if self.Nrep == 1:
                self._run_2D()
            else:
                if how == "sequential":
                    for i in range(self.Nrep):
                        self._run_2D(seed=i)
                elif how == "async":
                    with mp.Pool(mp.cpu_count()) as pool:
                        res = pool.starmap_async(
                            self._run_2D,
                            [(i,) for i in range(self.Nrep)],
                        ).get()
                        self._statsMultiple = {i: r for i, r in enumerate(res)}
                self._keys = [
                    "T_nuc_min",
                    "T_nuc_kin",
                    "T_nuc_mean",
                    "T_nuc_max",
                    "t_nuc",
                    "t_sol",
                    "t_fr",
                ]
        # update simulation status
        self._simulationStatus = 1

    # -------------------------------------------------------------------------- #
    # 0D implementation of the freezing model
    # -------------------------------------------------------------------------- #

    # homogeneous model
    def _run_0D(self, seed=0):
        # write into local vars for convenience: geometry
        A = self.const["A"]
        V = self.const["V"]

        # density
        rho_l = self.const["rho_l"]

        # mass of individual components
        mass = self.const["mass"]
        mass_water = self.const["mass_water"]
        mass_solute = self.const["mass_solute"]

        # solute mass fraction
        w_s = mass_solute / mass

        # specific heat capacities
        cp_w = self.const["cp_w"]
        cp_i = self.const["cp_i"]
        cp_s = self.const["cp_s"]
        cp_solution = self.const["cp_solution"]

        # soluton concentration and freezing-point depression
        solid_fraction = self.const["solid_fraction"]
        T_m = self.const["T_eq"] + 273.15
        k_f = self.const["k_f"]
        M_s = self.const["M_s"]
        depression = self.const["depression"]
        T_eq_l = T_m - depression

        # nucleation kinetics
        # kb = self.const["kb"]
        a = self.const["a"]
        b = self.const["b"]
        c = self.const["c"]

        # pre-exponential nucleation parameter
        np.random.seed(2024)
        xi_v = norm.ppf(np.random.rand())
        kb = 10 ** (-(a + xi_v * c))

        # latent heat of fusion
        Dh = self.const["Dh"]

        # heat transfer coefficient shelf
        K_shelf = self.k["s0"]

        # time step
        dt = 0.1  # [s]

        # initial temperature
        T_0 = self.opcond.cooling["start"] + 273.15
        T_k = T_0  # [K] temperature profile at t = T_old
        T_new = T_k

        # shelf temperature and coooling rate
        T_shelf_cool = self.opcond.tempProfile(dt) + 273.15

        # preallocation of results array
        T_product_cooling = np.zeros_like(T_shelf_cool)
        time_cooling = np.zeros_like(T_shelf_cool)

        # initial temperature
        T_old = T_0

        # random numbers
        np.random.seed(seed)
        F_rand = np.random.random()

        # preallocate the expected number of nuclei
        E_t, Nt_cool_end, Nt_sol_end = 0, None, None

        # solving for the cooling stage
        for i, T_shelf in enumerate(T_shelf_cool):
            # computing new temperature
            T_new = T_old + dt * (A * K_shelf * (T_shelf - T_old)) / (
                cp_solution * mass
            )
            # saving results
            T_product_cooling[i] = T_new
            time_cooling[i] = dt * i
            # update temperature
            T_old = T_new
            # check if nucleation occurs stochastically
            if T_eq_l > T_new:
                J = kb * (T_eq_l - T_new) ** b
            else:
                J = 0
            # frequency of nucleation events
            K_v = J * V
            # expected number of nuclei
            if K_v >= 0:
                E_t += K_v * dt
            # CDF of probability
            F_nuc = 1 - np.exp(-E_t)
            if self.opcond.cnTemp == None:
                # Monte-Carlo approach
                if F_nuc > F_rand:
                    Nt_cool_end = i
                    break
            else:
                # controlled nucleation
                if T_new <= (self.opcond.cnTemp + 273.15):
                    Nt_cool_end = i
                    break

        # check if nucleation occured
        if Nt_cool_end is None:
            raise ValueError("Nucleation did not occur, prolong process time.")

        # final solution arrays for cooling stage
        T_product_cooling = T_product_cooling[0:Nt_cool_end]
        T_shelf_cooling = T_shelf_cool[0:Nt_cool_end]
        time_cooling = time_cooling[0:Nt_cool_end]

        # mass fraction of ice over time
        w_i_cooling = np.zeros(len(T_shelf_cooling))

        # 2. Nucleation phase
        t_nuc = Nt_cool_end * dt
        T_nuc = T_new

        # save nucleation results
        self._stats["T_nuc"] = T_nuc - 273.15
        self._stats["t_nuc"] = t_nuc / 60
        self._stats["t_sol"] = None
        self._stats["t_fr"] = None

        # solving the quadratic equation
        B = -T_m - T_nuc - Dh * mass_water / (cp_solution * mass)
        C = (
            Dh * mass_water * T_m / (cp_solution * mass)
            - mass_solute * (k_f / M_s) * Dh / (cp_solution * mass)
            + T_m * T_nuc
        )

        # equilibrium temperatures determined from quadratic equation (A = 1)
        T_eq = 0.5 * (-B - np.sqrt(B**2 - 4 * C))

        # computed mass of ice from the quadratic equation
        m_i_eq = mass_water - mass_solute * (k_f / M_s) / (T_m - T_eq)

        # solution of ice nucleation stage
        w_i_nucl = m_i_eq / mass

        # preallocation
        T_product_solid = np.zeros(len(T_shelf_cool[Nt_cool_end:]))
        w_i_solid = np.zeros(len(T_shelf_cool[Nt_cool_end:]))
        time_solid = np.zeros(len(T_shelf_cool[Nt_cool_end:]))

        # saving results of the ice nucleation stage
        w_i_solid[0] = w_i_nucl
        T_product_solid[0] = T_eq

        # shelf temperatures
        T_old = T_eq
        w_i_new = w_i_nucl

        # solidification stage
        for i, T_shelf in enumerate(T_shelf_cool[Nt_cool_end:]):
            # thermal properties
            cp = cp_s * w_s + cp_i * w_i_new + cp_w * (1 - w_s - w_i_new)  # [J/kgK]
            # temperature computation
            T_new = T_old + (
                dt
                * (A * K_shelf * (T_shelf - T_old))
                * (
                    cp * rho_l * V
                    + (Dh * k_f * mass_solute / M_s) * (1 / (T_m - T_old) ** 2)
                )
                ** -1
            )
            # update temperatures
            T_old = T_new
            # save results
            T_product_solid[i] = T_new
            time_solid[i] = t_nuc + dt * i
            # new ice mass fraction
            w_i_new = (mass_water - (k_f * mass_solute / M_s) / (T_m - T_new)) / mass
            # solidification sigma
            sigma_new = w_i_new * mass / (mass - mass_solute)
            # save solidification time
            if (self._stats["t_sol"] is None) & (sigma_new >= 0.9):
                Nt_sol_end = i
                self._stats["t_sol"] = dt * i / 60
                self._stats["t_fr"] = (t_nuc + dt * i) / 60

        # check if solidification is completed
        if Nt_sol_end is None:
            raise ValueError("Solidification is not completed, prolong process time.")

        # ice mass fraction evolution
        w_i_solid = (
            mass_water - (k_f * mass_solute / M_s) / (T_m - T_product_solid)
        ) / mass

        # time
        time_results = np.concatenate((time_cooling, time_solid), axis=0) / 3600

        # temperatures
        temp_results = (
            np.concatenate((T_product_cooling, T_product_solid), axis=0) - 273.15
        )
        shelf_results = T_shelf_cool - 273.15

        # ice mass fraction
        ice_results = np.concatenate((w_i_cooling, w_i_solid), axis=0)

        # save results
        self._domain = 0
        self._prop = (T_eq_l, w_s)
        self._time = time_results
        self._shelfTemp = shelf_results
        self._temp = temp_results
        self._iceMassFraction = ice_results

        return [
            self._stats["T_nuc"],
            self._stats["t_nuc"],
            self._stats["t_sol"],
            self._stats["t_fr"],
        ]

    # -------------------------------------------------------------------------- #
    # 1D implementation of the freezing model
    # -------------------------------------------------------------------------- #

    # spatial model in 1 dimension
    def _run_1D(self, seed=0):
        # write into local vars for convenience: geometry
        height = self.const["height"]
        A = self.const["A"]
        V = self.const["V"]

        # density
        rho_l = self.const["rho_l"]

        # mass of individual components
        mass = self.const["mass"]
        mass_water = self.const["mass_water"]
        mass_solute = self.const["mass_solute"]
        w_s = mass_solute / mass

        # heat conductivities
        lambda_w = self.const["lambda_w"]
        lambda_i = self.const["lambda_i"]
        lambda_s = self.const["lambda_s"]

        # specific heat capacities
        cp_w = self.const["cp_w"]
        cp_i = self.const["cp_i"]
        cp_s = self.const["cp_s"]
        cp_solution = self.const["cp_solution"]

        # soluton concentration and freezing-point depression
        solid_fraction = self.const["solid_fraction"]
        T_m = self.const["T_eq"] + 273.15
        k_f = self.const["k_f"]
        M_s = self.const["M_s"]
        depression = self.const["depression"]
        T_eq_l = T_m - depression

        # nucleation kinetics
        # kb = self.const["kb"]
        a = self.const["a"]
        b = self.const["b"]
        c = self.const["c"]

        # pre-exponential nucleation parameter
        np.random.seed(2024)
        xi_v = norm.ppf(np.random.rand())
        kb = 10 ** (-(a + xi_v * c))

        # general
        k_B = self.const["k_B"]

        # latent heat of fusion
        Dh = self.const["Dh"]

        # heat transfer coefficient shelf
        K_shelf = self.k["s0"]

        # additional VISF parameters
        if self.const["configuration"] == "VISF":
            # vacuum pressure
            p_vac = self.const["p_vac"]
            # evaporation coefficient
            kappa = self.const["kappa"]
            # latent heat of vaporization for water [J/kg]
            dHe = self.const["Dh_evaporation"]
            # mass of water molecule [kg]
            m_water = self.const["m_water"]
            # time for vacuum start [h]
            t_vac_start = self.const["t_vac_start"]
            # duration of the vacuum [h]
            t_vac_duration = self.const["t_vac_duration"]

        # [W/mK] effective heat conductivity
        lambda_eff = solid_fraction * lambda_s + (1 - solid_fraction) * lambda_w
        # heat diffusivity
        diffusivity_cooling = lambda_eff / (cp_solution * rho_l)

        Nz = 30  # [/] number of spatial grid points
        dz = height / Nz  # [m] size of spatial step
        z = np.linspace(0, height, Nz)  # [m] z-position array

        # limiting value of alpha for computing dt from CFL condition
        alpha_max = lambda_i / (cp_i * rho_l)

        # temporal discretization
        dt = (
            0.4 * dz**2 / alpha_max
        )  # [s] size of time step (determined by the CFL cond.)

        Nt_exp = (
            int(np.ceil(self.opcond.t_tot / dt)) + 1
        )  # the total number of timesteps

        # Preallocation of arrays for saving results from the cooling stage.
        N_save_cool = 10000
        cooling_temp_results = np.zeros((N_save_cool, Nz))
        cooling_ice_results = np.zeros((N_save_cool, Nz))
        cooling_shelf_results = np.zeros(N_save_cool)
        cooling_time_results = np.zeros(N_save_cool)
        i_save = 0

        # initial temperature
        T_0 = self.opcond.cooling["start"] + 273.15
        T_k = np.zeros(Nz) + T_0  # [K] temperature profile at t = T_old

        # shelf temperature and coooling rate
        T_shelf_cool = self.opcond.tempProfile(dt) + 273.15

        # random numbers
        np.random.seed(seed)
        F_rand = np.random.random()

        # preallocate the expected number of nuclei
        E_t, Nt_cool_end, Nt_sol_end = 0, None, None

        for i, T_shelf in enumerate(T_shelf_cool):
            # overall heat flux
            q_overall = K_shelf * (T_shelf - T_k[0])
            # boundary condition on the bottom of the vial: ghost grid point
            T_bottom_BC = T_k[0] + q_overall * dz / lambda_eff

            if self.const["configuration"] == "VISF":
                # temperature at the top of the liquid
                T_l = T_k[-1]
                # vapour temperature
                T_v = T_k[-1]
                # vapour pressure
                p_vap = Utils.vapour_pressure_liquid(T_l)

                # evaporative heat flux (only in the span of vacuum, hold_2)
                if (dt * i > t_vac_start * 3600) & (
                    dt * i < (t_vac_start + t_vac_duration) * 3600
                ):
                    # vapour flux
                    N_w = Utils.vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v)
                    # evaporative heat flux
                    q_e = -N_w * dHe
                else:
                    q_e = 0

            else:
                q_e = 0

            # boundary condition in the top of the vial
            T_top_BC = T_k[Nz - 1] + q_e * dz / lambda_eff

            # updating the temperature at the bottom
            T_bottom = T_k[0] + (diffusivity_cooling * dt / dz**2) * (
                T_k[1] - 2 * T_k[0] + T_bottom_BC
            )

            # updating the bulk temperatures
            T_center = T_k[1 : Nz - 1] + (diffusivity_cooling * dt / dz**2) * (
                T_k[2:Nz] - 2 * T_k[1 : Nz - 1] + T_k[0 : Nz - 2]
            )

            # Top point.
            # updating the top temperature
            T_top = T_k[Nz - 1] + (diffusivity_cooling * dt / dz**2) * (
                T_top_BC - 2 * T_k[Nz - 1] + T_k[Nz - 2]
            )

            # update temperatures
            T_k = np.concatenate((T_bottom, T_center, T_top), axis=None)

            # nucleation rate (global aproach)
            superCooledMask = T_k < T_eq_l
            J_z = np.zeros(z.shape)
            J_z[superCooledMask] = kb * (T_eq_l - T_k[superCooledMask]) ** b

            # compute the frequency of nucleation
            K_v = A * simps(J_z, z)

            # check the range of LCS = LocalSupercooling
            LCS_i = T_k < T_eq_l
            # the not so supercool part
            LCS_i_r = ~LCS_i

            # sparsely saving results
            if np.mod(i, np.ceil(Nt_exp / N_save_cool)) == 0:
                # saving temperature results
                cooling_temp_results[i_save, :] = T_k - 273.15
                # saving shelf temperature results
                cooling_shelf_results[i_save] = T_shelf - 273.15
                # ice
                cooling_ice_results[i_save] = 0 * T_k
                # saving sparse time array
                cooling_time_results[i_save] = dt * i
                # update saving index
                i_save = i_save + 1

            # expected number of nuclei
            E_t += K_v * dt
            # CDF of probability
            F_nuc = 1 - np.exp(-E_t)
            if self.opcond.cnTemp == None:
                # Monte-Carlo approach
                if F_nuc > F_rand:
                    Nt_cool_end = i
                    t_nuc = dt * i
                    T_nuc = T_k
                    break
            else:
                # controlled nucleation
                if T_k.any() <= (self.opcond.cnTemp + 273.15):
                    Nt_cool_end = i
                    t_nuc = dt * i
                    T_nuc = T_k
                    break

        # check if nucleation occured
        if Nt_cool_end is None:
            raise ValueError("Nucleation did not occur, prolong process time.")

        # kinetic mean nucleation temperature
        if K_v > 0:
            T_nuc_kin = (A / K_v) * simps(T_nuc * J_z, z)
        else:
            T_nuc_kin = 273.15

        # save nucleation temperatures
        self._stats = {
            "T_nuc_min": T_nuc.min() - 273.15,
            "T_nuc_kin": T_nuc_kin - 273.15,
            "T_nuc_mean": np.mean(T_nuc) - 273.15,
            "T_nuc_max": T_nuc.max() - 273.15,
        }

        # save nucleation time
        self._stats["t_nuc"] = t_nuc / 60
        self._stats["t_sol"] = None
        self._stats["t_fr"] = None

        " 2. Ice nucleation stage. "

        # solving the quadratic equation (vectorized (A = 1)
        B = -T_m - T_nuc - Dh * mass_water / (cp_solution * mass)
        C = (
            Dh * mass_water * T_m / (cp_solution * mass)
            - mass_solute * (k_f / M_s) * Dh / (cp_solution * mass)
            + T_m * T_nuc
        )

        # lcoal supercooling
        LCS_i = T_k < T_eq_l
        # the not so supercool part
        LCS_i_r = ~LCS_i

        # equilibrium temperatures determined from quadratic equation (A = 1)
        T_eq_sol_1 = 0.5 * (-B - np.sqrt(B**2 - 4 * C))

        # temperaure field after nucleation
        T_after_nuc = T_k * LCS_i_r + T_eq_sol_1 * LCS_i

        # computed mass of ice from the quadratic equation (T_k used just to get
        # the same dimensions)
        m_i_nucl_sol_1 = mass_water - mass_solute * (k_f / M_s) / (T_m - T_eq_sol_1)
        m_i_eq = 0 * T_k * LCS_i_r + m_i_nucl_sol_1 * LCS_i

        # save temperature profile after nucleation
        cooling_temp_results[i_save, :] = T_after_nuc - 273.15
        # saving shelf temperature results
        cooling_shelf_results[i_save] = T_shelf - 273.15
        # ice
        cooling_ice_results[i_save] = m_i_eq / (mass_water + mass_solute)
        # saving sparse time array
        cooling_time_results[i_save] = dt * i

        i_end = i
        i_save_end = i_save

        " 3. Solidification stage. "

        N_save_solid = 10000

        # arrays of temperature
        T_k = T_after_nuc  # [K] temperature profile at t = T_old

        # array of ice mass fractions
        w_i_k = m_i_eq / (mass_water + mass_solute)
        Nt_solid = Nt_exp - i_end

        # Preallocation of arrays for saving results from the cooling stage.
        solid_temp_results = np.zeros((N_save_solid, Nz))
        solid_ice_results = np.zeros((N_save_solid, Nz))
        solid_shelf_results = np.zeros(N_save_solid)
        solid_time_results = np.zeros(N_save_solid)
        i_save = 0

        # FDM scheme for solving the heat conduction equation
        for i, T_shelf in enumerate(T_shelf_cool[i_end:]):
            # local supercooling
            LCS_i = T_k < T_eq_l
            # the not so supercool part
            LCS_i_r = ~LCS_i

            " Thermal properties. "
            # heat capacity
            cp_solution = (
                cp_s * solid_fraction
                + cp_i * w_i_k
                + cp_w * (1 - solid_fraction - w_i_k)
            )
            # heat conductivity
            lambda_eff = lambda_i * w_i_k + lambda_w * (1 - w_i_k)
            # nonlinear capacitance term
            beta = Dh * k_f * mass_solute / (M_s * rho_l * V * cp_solution)
            # total nonlinear capacitance term
            BETA = np.ones(Nz) * LCS_i_r + (1 + beta / (T_k - T_m) ** 2) * LCS_i

            " Boundary condition at the bottom, z = 0."
            # overall heat flux
            q_overall = K_shelf * (T_shelf - T_k[0])
            # boundary condition on the bottom of the vial: ghost grid point
            T_bottom_BC = T_k[0] + q_overall * dz / lambda_eff[0]

            if self.const["configuration"] == "VISF":
                # temperature at the top of the liquid
                T_l = T_k[-1]
                # vapour temperature
                T_v = T_k[-1]
                # vapour pressure
                p_vap = Utils.vapour_pressure_solid(T_l)

                # evaporative heat flux (only in the span of vacuum, hold_2)
                if (t_nuc + dt * i > t_vac_start * 3600) & (
                    t_nuc + dt * i < (t_vac_start + t_vac_duration) * 3600
                ):
                    # vapour flux
                    N_w = Utils.vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v)
                    # evaporative heat flux
                    q_e = -N_w * dHe
                else:
                    q_e = 0

            else:
                q_e = 0

            " Boundary condition at the top, z = height."
            # boundary condition in the top of the vial
            T_top_BC = T_k[Nz - 1] + q_e * dz / lambda_eff[Nz - 1]

            " Bottom point. "
            # update the point at the bottom
            T_bottom = T_k[0] + (dt / (cp_solution[0] * rho_l)) * (
                (lambda_eff[1] - lambda_eff[0]) * (T_k[1] - T_bottom_BC) / (4 * dz**2)
                + lambda_eff[0] * (T_k[1] - 2 * T_k[0] + T_bottom_BC) / dz**2
            ) * BETA[0] ** (-1)

            " Center points. "
            # update the bulk points
            T_center = T_k[1 : Nz - 1] + (dt / (cp_solution[1 : Nz - 1] * rho_l)) * (
                (lambda_eff[2:Nz] - lambda_eff[0 : Nz - 2])
                * (T_k[2:Nz] - T_k[0 : Nz - 2])
                / (4 * dz**2)
                + lambda_eff[1 : Nz - 1]
                * (T_k[2:Nz] - 2 * T_k[1 : Nz - 1] + T_k[0 : Nz - 2])
                / dz**2
            ) * BETA[1 : Nz - 1] ** (-1)

            " Top point. "
            # update the temperature at the top
            T_top = T_k[Nz - 1] + (dt / (cp_solution[Nz - 1] * rho_l)) * (
                (lambda_eff[Nz - 1] - lambda_eff[Nz - 2])
                * (T_top_BC - T_k[Nz - 2])
                / (4 * dz**2)
                + lambda_eff[Nz - 1]
                * (T_top_BC - 2 * T_k[Nz - 1] + T_k[Nz - 2])
                / dz**2
            ) * BETA[Nz - 1] ** (-1)

            " Processing the results I: temperatures. "
            # update temperatures
            T_k = np.concatenate((T_bottom, T_center, T_top), axis=None)

            " Check for local supercooling. "
            # check the range of LCS
            LCS_i = np.where(T_k < T_eq_l, 1, 0)
            # the not so supercool region
            LCS_i_r = np.where(LCS_i, 0, 1)

            " Processing the results II: ice mass fractions. "
            # ice mass distribution
            m_ice = (
                np.zeros(Nz) * LCS_i_r
                + (mass_water - mass_solute * (k_f / M_s) / (T_m - T_k)) * LCS_i
            )
            # ice mass fraction distribution
            w_i_k = m_ice / mass

            " Saving results of the cooling stage. "
            # sparse saving

            if np.mod(i, np.ceil(Nt_solid / N_save_solid)) == 0:
                # temperature results
                solid_temp_results[i_save, :] = T_k - 273.15
                # ice mass fraction results
                solid_ice_results[i_save, :] = w_i_k
                # shelf temperature results
                solid_shelf_results[i_save] = T_shelf - 273.15
                # saving sparse time array
                solid_time_results[i_save] = t_nuc + dt * i
                # updating saving index
                i_save = i_save + 1

            # solidification sigma
            sigma_new = (1 / z[-1]) * simps(m_ice, z) / (mass - mass_solute)
            # save solidification time
            if (self._stats["t_sol"] is None) & (sigma_new >= 0.9):
                Nt_sol_end = i
                self._stats["t_sol"] = dt * i / 60
                self._stats["t_fr"] = (t_nuc + dt * i) / 60

        # check if solidification is completed
        if Nt_sol_end is None:
            raise ValueError("Solidification is not completed, prolong process time.")

        # temperature distributions in the product
        temp_results = np.concatenate(
            (
                cooling_temp_results[: i_save_end + 1, :],
                solid_temp_results[: i_save - 1, :],
            ),
            axis=0,
        )

        # ice mass fraction results
        ice_results = np.concatenate(
            (
                cooling_ice_results[: i_save_end + 1, :],
                solid_ice_results[: i_save - 1, :],
            ),
            axis=0,
        )

        # total time array
        time_results = (
            np.concatenate(
                (
                    cooling_time_results[: i_save_end + 1],
                    solid_time_results[: i_save - 1],
                ),
                axis=0,
            )
            / 3600
        )

        # total shelf temperature array
        shelf_results = np.concatenate(
            (
                cooling_shelf_results[: i_save_end + 1],
                solid_shelf_results[: i_save - 1],
            ),
            axis=0,
        )

        # save results
        self._domain = z
        self._nuc = (t_nuc, T_nuc)
        self._prop = (T_eq_l, w_s)
        self._time = time_results
        self._shelfTemp = shelf_results
        self._temp = temp_results
        self._iceMassFraction = ice_results

        return [
            self._stats["T_nuc_min"],
            self._stats["T_nuc_kin"],
            self._stats["T_nuc_mean"],
            self._stats["T_nuc_max"],
            self._stats["t_nuc"],
            self._stats["t_sol"],
            self._stats["t_fr"],
        ]

    # -------------------------------------------------------------------------- #
    # 2D implementation of the freezing model
    # -------------------------------------------------------------------------- #

    # spatial model in 2 dimensions
    def _run_2D(self, seed=0):
        # write into local vars for convenience: geometry
        height = self.const["height"]
        diameter = self.const["diameter"]
        V = self.const["V"]
        radius = diameter / 2

        # density
        rho_l = self.const["rho_l"]

        # mass of individual components
        mass = self.const["mass"]
        mass_water = self.const["mass_water"]
        mass_solute = self.const["mass_solute"]
        w_s = mass_solute / mass

        # heat conductivities
        lambda_w = self.const["lambda_w"]
        lambda_i = self.const["lambda_i"]
        lambda_s = self.const["lambda_s"]

        # specific heat capacities
        cp_w = self.const["cp_w"]
        cp_i = self.const["cp_i"]
        cp_s = self.const["cp_s"]
        cp_solution = self.const["cp_solution"]

        # soluton concentration and freezing-point depression
        solid_fraction = self.const["solid_fraction"]
        T_m = self.const["T_eq"] + 273.15
        k_f = self.const["k_f"]
        M_s = self.const["M_s"]
        depression = self.const["depression"]
        T_eq_l = T_m - depression

        # nucleation kinetics
        # kb = self.const["kb"]
        a = self.const["a"]
        b = self.const["b"]
        c = self.const["c"]

        # pre-exponential nucleation parameter
        np.random.seed(2024)
        xi_v = norm.ppf(np.random.rand())
        kb = 10 ** (-(a + xi_v * c))

        # general
        k_B = self.const["k_B"]

        # latent heat of fusion
        Dh = self.const["Dh"]

        # shelf heat transfer coefficient
        K_shelf = self.k["s0"]

        # additional VISF parameters
        if self.const["configuration"] == "VISF":
            # vacuum pressure
            p_vac = self.const["p_vac"]
            # evaporation coefficient
            kappa = self.const["kappa"]
            # latent heat of vaporization for water [J/kg]
            dHe = self.const["Dh_evaporation"]
            # mass of water molecule [kg]
            m_water = self.const["m_water"]
            # time for vacuum start [h]
            t_vac_start = self.const["t_vac_start"]
            # duration of the vacuum [h]
            t_vac_duration = self.const["t_vac_duration"]

        # additional jacket-ramped freezing parameters
        if self.const["configuration"] == "jacket":
            # vacuum pressure
            air_gap = self.const["air_gap"]
            # evaporation coefficient
            lambda_air = self.const["lambda_air"]
            # compute the wall heat transfer coefficient
            K_wall = (1 / K_shelf + air_gap / lambda_air) ** (-1)  # [W/m2K]

        # [W/mK] effective heat conductivity
        k_eff = solid_fraction * lambda_s + (1 - solid_fraction) * lambda_w

        # heat diffusivity
        alpha = k_eff / (cp_solution * rho_l)

        Nz = 30  # [/] number of spatial grid points
        dz = height / Nz  # [m] size of spatial step
        z = np.linspace(0, height, Nz)  # [m] z-position array

        # spatal discretization in the radial direction
        Nr = 15  # [/] number of spatial grid
        dr = radius / Nr  # [m] size of spatial step
        r = np.linspace(0, radius, Nr)  # [m] z-position array

        # limiting value of alpha for computing dt from CFL condition
        alpha_max = lambda_i / (cp_i * rho_l)

        # temporal discretization
        dt = (
            (0.4 / alpha_max) * (dz**2 * dr**2) / (dr**2 + dz**2)
        )  # [s] size of time step

        Nt_exp = (
            int(np.ceil(self.opcond.t_tot / dt)) + 1
        )  # the total number of timesteps

        # Preallocation of arrays for saving results from the cooling stage.
        N_save_cool = 10000
        cooling_temp_results = np.zeros((N_save_cool, Nz, Nr))
        cooling_ice_results = np.zeros((N_save_cool, Nz, Nr))
        cooling_shelf_results = np.zeros(N_save_cool)
        cooling_time_results = np.zeros(N_save_cool)
        i_save = 0

        # initial temperature
        T_0 = self.opcond.cooling["start"] + 273.15
        T_k = np.zeros((Nz, Nr)) + T_0  # [K] temperature profile at t = T_old
        T_new = T_k

        # shelf temperature and coooling rate
        T_shelf_cool = self.opcond.tempProfile(dt) + 273.15

        # random numbers
        np.random.seed(seed)
        F_rand = np.random.random()

        # preallocate the expected number of nuclei
        E_t, Nt_cool_end, Nt_sol_end = 0, None, None

        # finite difference method for solving the heat conduction equation
        for i, T_shelf in enumerate(T_shelf_cool):
            "Boundary condition at the bottom, z = 0."
            # overall heat flux
            q_overall = K_shelf * (T_shelf - T_k[0, :])
            # boundary condition on the bottom of the vial: ghost grid point
            T_bottom = T_k[0, :] + q_overall * dz / k_eff

            if self.const["configuration"] == "VISF":
                # temperature at the top of the liquid
                T_l = T_k[-1, :]
                # vapour temperature
                T_v = T_k[-1, :]
                # vapour pressure
                p_vap = Utils.vapour_pressure_solid(T_l)

                # evaporative heat flux (only in the span of vacuum, hold_2)
                if (dt * i > t_vac_start * 3600) & (
                    dt * i < (t_vac_start + t_vac_duration) * 3600
                ):
                    # vapour flux
                    N_w = Utils.vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v)
                    # evaporative heat flux
                    q_e = -N_w * dHe
                else:
                    q_e = 0

            else:
                q_e = 0

            # boundary condition in the top of the vial
            T_top = T_k[Nz - 1, :] + q_e * dz / k_eff

            " Boundary conditions at the edge, r = R. "

            if self.const["configuration"] == "jacket":
                # jacket heat flux
                q_jacket = K_wall * (T_shelf - T_k[:, Nr - 1])

            else:
                # no cooling from jacket
                q_jacket = 0

            # boundary condition in the top of the vial
            T_edge = T_k[:, Nr - 1] + q_jacket * dz / k_eff

            " Boundary condition in the centre, r = 0. "
            # axial symmetry
            T_center = T_k[:, 0]

            " Bottom boundary points: z = 0. "
            # bottom-centre: r = 0
            T_new[0, 0] = T_k[0, 0] + alpha * dt * (
                2 * (T_k[0, 1] - 2 * T_k[0, 0] + T_center[0]) / dr**2
                + (T_k[1, 0] - 2 * T_k[0, 0] + T_bottom[0]) / dz**2
            )

            # bottom corner (left or right - symmetry): r = R
            T_new[0, Nr - 1] = T_k[0, Nr - 1] + alpha * dt * (
                (1 / r[Nr - 1]) * (T_edge[0] - T_k[0, Nr - 2]) / (2 * dr)
                + (T_edge[0] - 2 * T_k[0, Nr - 1] + T_k[0, Nr - 2]) / dr**2
                + (T_k[1, Nr - 1] - 2 * T_k[0, Nr - 1] + T_bottom[Nr - 1]) / dz**2
            )

            # the rest of the bottom: 0 < r < R
            T_new[0, 1 : Nr - 1] = T_k[0, 1 : Nr - 1] + alpha * dt * (
                (1 / r[1 : Nr - 1]) * (T_k[0, 2:Nr] - T_k[0, 0 : Nr - 2]) / (2 * dr)
                + (T_k[0, 2:Nr] - 2 * T_k[0, 1 : Nr - 1] + T_k[0, 0 : Nr - 2]) / dr**2
                + (T_k[1, 1 : Nr - 1] - 2 * T_k[0, 1 : Nr - 1] + T_bottom[1 : Nr - 1])
                / dz**2
            )

            " Top boundary points: z = H. "
            # top centre: r = 0
            T_new[Nz - 1, 0] = T_k[Nz - 1, 0] + alpha * dt * (
                2 * (T_k[Nz - 1, 1] - 2 * T_k[Nz - 1, 0] + T_center[-1]) / dr**2
                + (T_top[0] - 2 * T_k[Nz - 1, 0] + T_k[Nz - 2, 0]) / dz**2
            )

            # top corner (left and right - symmetry): r = R
            T_new[Nz - 1, Nr - 1] = T_k[Nz - 1, Nr - 1] + alpha * dt * (
                (1 / r[Nr - 1]) * (T_edge[Nz - 1] - T_k[Nz - 1, Nr - 2]) / (2 * dr)
                + (T_edge[Nz - 1] - 2 * T_k[Nz - 1, Nr - 1] + T_k[Nz - 1, Nr - 2])
                / dr**2
                + (T_top[Nr - 1] - 2 * T_k[Nz - 1, Nr - 1] + T_k[Nz - 2, Nr - 1])
                / dz**2
            )

            # the rest of the vial's top: 0 < r < R
            T_new[Nz - 1, 1 : Nr - 1] = T_k[Nz - 1, 1 : Nr - 1] + alpha * dt * (
                (1 / r[1 : Nr - 1])
                * (T_k[Nz - 1, 2:Nr] - T_k[Nz - 1, 0 : Nr - 2])
                / (2 * dr)
                + (
                    T_k[Nz - 1, 2:Nr]
                    - 2 * T_k[Nz - 1, 1 : Nr - 1]
                    + T_k[Nz - 1, 0 : Nr - 2]
                )
                / dr**2
                + (
                    T_top[1 : Nr - 1]
                    - 2 * T_k[Nz - 1, 1 : Nr - 1]
                    + T_k[Nz - 2, 1 : Nr - 1]
                )
                / dz**2
            )

            " Edge boundary points: r = R. "
            # curved cylinder surface in contact with air: 0 < z < H
            T_new[1 : Nz - 1, Nr - 1] = T_k[1 : Nz - 1, Nr - 1] + alpha * dt * (
                (1 / r[Nr - 1])
                * (T_edge[1 : Nz - 1] - T_k[1 : Nz - 1, Nr - 2])
                / (2 * dr)
                + (
                    T_edge[1 : Nz - 1]
                    - 2 * T_k[1 : Nz - 1, Nr - 1]
                    + T_k[1 : Nz - 1, Nr - 2]
                )
                / dr**2
                + (
                    T_k[2:Nz, Nr - 1]
                    - 2 * T_k[1 : Nz - 1, Nr - 1]
                    + T_k[0 : Nz - 2, Nr - 1]
                )
                / dz**2
            )

            " Center boundary points: r = 0. "
            # symmetry boundary condition in the centre of the cylinder: 0 < z < H
            T_new[1 : Nz - 1, 0] = T_k[1 : Nz - 1, 0] + alpha * dt * (
                2
                * (T_k[1 : Nz - 1, 1] - 2 * T_k[1 : Nz - 1, 0] + T_center[1 : Nz - 1])
                / dr**2
                + (T_k[2:Nz, 0] - 2 * T_k[1 : Nz - 1, 0] + T_k[0 : Nz - 2, 0]) / dz**2
            )

            " Bulk of the vial points: : 0 < z < H & 0 < r < R. "
            # all the other points in the domain
            T_new[1 : Nz - 1, 1 : Nr - 1] = T_k[1 : Nz - 1, 1 : Nr - 1] + alpha * dt * (
                (1 / r[1 : Nr - 1])
                * (T_k[1 : Nz - 1, 2:Nr] - T_k[1 : Nz - 1, 0 : Nr - 2])
                / (2 * dr)
                + (
                    T_k[1 : Nz - 1, 2:Nr]
                    - 2 * T_k[1 : Nz - 1, 1 : Nr - 1]
                    + T_k[1 : Nz - 1, 0 : Nr - 2]
                )
                / dr**2
                + (
                    T_k[2:Nz, 1 : Nr - 1]
                    - 2 * T_k[1 : Nz - 1, 1 : Nr - 1]
                    + T_k[0 : Nz - 2, 1 : Nr - 1]
                )
                / dz**2
            )

            " Processing the results of the time iteration. "
            # update temperatures
            T_k = T_new

            # Stochastic nucleation for the entire vial.
            # nucleation rate (global aproach)
            superCooledMask = T_k < T_eq_l
            J_z_r = np.zeros((Nz, Nr))
            J_z_r[superCooledMask] = kb * (T_eq_l - T_k[superCooledMask]) ** b

            # frequency of nucleation events itegrated over z
            K_r = 2 * np.pi * simps(r * J_z_r, r)
            # frequency of nucleation events itegrated over r
            K_v = simps(K_r, z)

            # Check for local supercooling.
            # check the range of LCS = LocalSupercooling
            LCS_i = T_k < T_eq_l
            # the not so supercool part
            LCS_i_r = ~LCS_i

            # Saving results of the cooling stage.
            # sparsely saving results
            if np.mod(i, np.ceil(Nt_exp / N_save_cool)) == 0:
                # saving temperature results
                cooling_temp_results[i_save, :, :] = T_k - 273.15
                # saving shelf temperature results
                cooling_shelf_results[i_save] = T_shelf - 273.15
                # ice
                cooling_ice_results[i_save, :, :] = 0 * T_k
                # saving sparse time array
                cooling_time_results[i_save] = dt * i
                # update saving index
                i_save = i_save + 1

            # expected number of nuclei
            E_t += K_v * dt
            # CDF of probability
            F_nuc = 1 - np.exp(-E_t)
            if self.opcond.cnTemp == None:
                # Monte-Carlo approach
                if F_nuc > F_rand:
                    Nt_cool_end = i
                    t_nuc = dt * i
                    T_nuc = T_k
                    break
            else:
                # controlled nucleation
                if T_k.any() <= (self.opcond.cnTemp + 273.15):
                    Nt_cool_end = i
                    t_nuc = dt * i
                    T_nuc = T_k
                    break

        # check if nucleation occured
        if Nt_cool_end is None:
            raise ValueError("Nucleation did not occur, prolong process time.")

        # compute the mean kinetic mean nucleation temperature
        f_z = 2 * np.pi * simps(r * T_nuc * J_z_r, r)
        if K_v > 0:
            T_nuc_kin = (1 / K_v) * simps(f_z, z)
        else:
            T_nuc_kin = 273.15

        # save nucleation temperatures
        self._stats = {
            "T_nuc_min": T_nuc.min() - 273.15,
            "T_nuc_kin": T_nuc_kin - 273.15,
            "T_nuc_mean": np.mean(T_nuc) - 273.15,
            "T_nuc_max": T_nuc.max() - 273.15,
        }

        # save nucleation time
        self._stats["t_nuc"] = t_nuc / 60
        self._stats["t_sol"] = None
        self._stats["t_fr"] = None

        # solving the quadratic equation (vectorized (A = 1)
        B = -T_m - T_nuc - Dh * mass_water / (cp_solution * mass)
        C = (
            Dh * mass_water * T_m / (cp_solution * mass)
            - mass_solute * (k_f / M_s) * Dh / (cp_solution * mass)
            + T_m * T_nuc
        )

        # lcoal supercooling
        LCS_i = T_k < T_eq_l
        # the not so supercool part
        LCS_i_r = ~LCS_i

        # equilibrium temperatures determined from quadratic equation (A = 1)
        T_eq_sol_1 = 0.5 * (-B - np.sqrt(B**2 - 4 * C))

        # temperaure field after nucleation
        T_after_nuc = T_k * LCS_i_r + T_eq_sol_1 * LCS_i

        # computed mass of ice from the quadratic equation (T_k used just to get
        # the same dimensions)
        m_i_nucl_sol_1 = mass_water - mass_solute * (k_f / M_s) / (T_m - T_eq_sol_1)
        m_i_eq = 0 * T_k * LCS_i_r + m_i_nucl_sol_1 * LCS_i

        cooling_temp_results[i_save, :, :] = T_after_nuc - 273.15
        # saving shelf temperature results
        cooling_shelf_results[i_save] = T_shelf - 273.15
        # ice
        cooling_ice_results[i_save] = m_i_eq / (mass_water + mass_solute)
        # saving sparse time array
        cooling_time_results[i_save] = dt * i
        i_end = i
        i_save_end = i_save
        N_save_solid = 10000

        # arrays of temperature
        T_k = T_after_nuc  # [K] temperature profile at t = T_old

        # array of ice mass fractions
        w_i_new = m_i_eq / (mass_water + mass_solute)

        Nt_solid = Nt_exp - i_end

        # Preallocation of arrays for saving results from the cooling stage.
        solid_temp_results = np.zeros((N_save_solid, Nz, Nr))
        solid_ice_results = np.zeros((N_save_solid, Nz, Nr))
        solid_shelf_results = np.zeros(N_save_solid)
        solid_time_results = np.zeros(N_save_solid)
        i_save = 0

        # FDM scheme for solving the heat conduction equation
        for i, T_shelf in enumerate(T_shelf_cool[i_end:]):
            # thermal properties
            # heat capacity
            cp_eff = (
                cp_s * solid_fraction
                + cp_i * w_i_new
                + cp_w * (1 - solid_fraction - w_i_new)
            )
            # heat conductivity
            k_eff = lambda_i * w_i_new + lambda_w * (1 - w_i_new)
            # heat diffusivity
            alpha = k_eff / (cp_eff * rho_l)
            # beta parameter beta = beta(t,z,r)
            beta = Dh * k_f * mass_solute / (M_s * rho_l * V * cp_eff)
            # the total nonlinear capacitance term
            BETA = np.ones((Nz, Nr)) * LCS_i_r + (1 + beta / (T_k - T_m) ** 2) * LCS_i

            # Boundary condition at the bottom, z = 0.
            # overall heat flux
            q_overall = K_shelf * (T_shelf - T_k[0, :])
            # q_overall = sigma_B*(T_shelf**4 - T_old[0,:]**4)
            # boundary condition on the bottom of the vial: ghost grid point
            T_bottom = T_k[0, :] + q_overall * dz / k_eff[0, :]

            if self.const["configuration"] == "VISF":
                # temperature at the top of the liquid
                T_l = T_k[-1, :]
                # vapour temperature
                T_v = T_k[-1, :]
                # vapour pressure
                p_vap = Utils.vapour_pressure_solid(T_l)

                # evaporative heat flux (only in the span of vacuum, hold_2)
                if (t_nuc + dt * i > t_vac_start * 3600) & (
                    t_nuc + dt * i < (t_vac_start + t_vac_duration) * 3600
                ):
                    # vapour flux
                    N_w = Utils.vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v)
                    # evaporative heat flux
                    q_e = -N_w * dHe
                else:
                    q_e = 0

            else:
                q_e = 0

            # boundary condition in the top of the vial
            T_top = T_k[Nz - 1, :] + q_e * dz / k_eff[Nz - 1, :]

            # boundary condition at the side of the vial
            if self.const["configuration"] == "jacket":
                # jacket heat flux
                q_jacket = K_wall * (T_shelf - T_k[:, Nr - 1])

            else:
                # no cooling from jacket
                q_jacket = 0

            # boundary condition in the top of the vial
            T_edge = T_k[:, Nr - 1] + q_jacket * dz / k_eff[:, Nr - 1]

            # Boundary condition in the centre, r = 0.
            # axial symmetry
            T_center = T_k[:, 0]

            # Bottom boundary points: z = 0.
            # bottom-centre: r = 0
            T_new[0, 0] = T_k[0, 0] + (dt / (cp_eff[0, 0] * rho_l)) * (
                2 * k_eff[0, 0] * (T_k[0, 1] - 2 * T_k[0, 0] + T_center[0]) / dr**2
                + (k_eff[0, 1] - k_eff[0, 0]) * (T_k[0, 1] - T_center[0]) / (4 * dr**2)
                + (k_eff[1, 0] - k_eff[0, 0]) * (T_k[1, 0] - T_bottom[0]) / (4 * dz**2)
                + k_eff[0, 0] * (T_k[1, 0] - 2 * T_k[0, 0] + T_bottom[0]) / dz**2
            ) * BETA[0, 0] ** (-1)

            # bottom corner (left or right - symmetry): r = R
            T_new[0, Nr - 1] = T_k[0, Nr - 1] + (dt / (cp_eff[0, Nr - 1] * rho_l)) * (
                (k_eff[0, Nr - 1] / r[Nr - 1]) * (T_edge[0] - T_k[0, Nr - 2]) / (2 * dr)
                + (k_eff[0, Nr - 1] - k_eff[0, Nr - 2])
                * (T_edge[0] - T_k[0, Nr - 1])
                / (4 * dr**2)
                + (k_eff[1, Nr - 1] - k_eff[0, Nr - 1])
                * (T_k[1, Nr - 1] - T_bottom[Nr - 1])
                / (4 * dz**2)
                + k_eff[0, Nr - 1]
                * (T_edge[0] - 2 * T_k[0, Nr - 1] + T_k[0, Nr - 2])
                / dr**2
                + k_eff[0, Nr - 1]
                * (T_k[1, Nr - 1] - 2 * T_k[0, Nr - 1] + T_bottom[Nr - 1])
                / dz**2
            ) * BETA[0, Nr - 1] ** (-1)

            # the rest of the bottom: 0 < r < R
            T_new[0, 1 : Nr - 1] = T_k[0, 1 : Nr - 1] + (
                dt / (cp_eff[0, 1 : Nr - 1] * rho_l)
            ) * (
                (k_eff[0, 1 : Nr - 1] / r[1 : Nr - 1])
                * (T_k[0, 2:Nr] - T_k[0, 0 : Nr - 2])
                / (2 * dr)
                + (k_eff[0, 2:Nr] - k_eff[0, 0 : Nr - 2])
                * (T_k[0, 2:Nr] - T_k[0, 0 : Nr - 2])
                / (4 * dr**2)
                + (k_eff[1, 1 : Nr - 1] - k_eff[0, 1 : Nr - 1])
                * (T_k[1, 1 : Nr - 1] - T_bottom[1 : Nr - 1])
                / (4 * dz**2)
                + k_eff[0, 1 : Nr - 1]
                * (T_k[0, 2:Nr] - 2 * T_k[0, 1 : Nr - 1] + T_k[0, 0 : Nr - 2])
                / dr**2
                + k_eff[0, 1 : Nr - 1]
                * (T_k[1, 1 : Nr - 1] - 2 * T_k[0, 1 : Nr - 1] + T_bottom[1 : Nr - 1])
                / dz**2
            ) * BETA[
                0, 1 : Nr - 1
            ] ** (
                -1
            )

            # Top boundary points: z = H.
            # top centre: r = 0
            T_new[Nz - 1, 0] = T_k[Nz - 1, 0] + (dt / (cp_eff[Nz - 1, 0] * rho_l)) * (
                2
                * k_eff[Nz - 1, 0]
                * (T_k[Nz - 1, 1] - 2 * T_k[Nz - 1, 0] + T_center[Nz - 1])
                / dr**2
                + (k_eff[Nz - 1, 1] - k_eff[Nz - 1, 0])
                * (T_k[Nz - 1, 1] - T_center[Nz - 1])
                / (4 * dr**2)
                + (k_eff[Nz - 1, 0] - k_eff[Nz - 2, 0])
                * (T_top[0] - T_k[Nz - 2, 0])
                / (4 * dz**2)
                + k_eff[Nz - 1, 0]
                * (T_top[0] - 2 * T_k[Nz - 1, 0] + T_k[Nz - 2, 0])
                / dz**2
            ) * BETA[Nz - 1, 0] ** (-1)

            # top edge corner (left or right - symmetry): r = R
            T_new[Nz - 1, Nr - 1] = T_k[Nz - 1, Nr - 1] + (
                dt / (cp_eff[Nz - 1, Nr - 1] * rho_l)
            ) * (
                (k_eff[Nz - 1, Nr - 1] / r[Nr - 1])
                * (T_edge[Nz - 1] - T_k[Nz - 1, Nr - 2])
                / (2 * dr)
                + (k_eff[Nz - 1, Nr - 1] - k_eff[Nz - 1, Nr - 2])
                * (T_edge[Nz - 1] - T_k[Nz - 1, Nr - 2])
                / (4 * dr**2)
                + (k_eff[Nz - 1, Nr - 1] - k_eff[Nz - 2, Nr - 1])
                * (T_top[Nr - 1] - T_k[Nz - 2, Nr - 1])
                / (4 * dz**2)
                + k_eff[Nz - 1, Nr - 1]
                * (T_edge[Nz - 1] - 2 * T_k[Nz - 1, Nr - 1] + T_k[Nz - 1, Nr - 2])
                / dr**2
                + k_eff[Nz - 1, Nr - 1]
                * (T_top[Nr - 1] - 2 * T_k[Nz - 1, Nr - 1] + T_k[Nz - 2, Nr - 1])
                / dz**2
            ) * BETA[
                Nz - 1, Nr - 1
            ] ** (
                -1
            )

            # the rest of the vial's top: 0 < r < R
            T_new[Nz - 1, 1 : Nr - 1] = T_k[Nz - 1, 1 : Nr - 1] + (
                dt / (cp_eff[Nz - 1, 1 : Nr - 1] * rho_l)
            ) * (
                (k_eff[Nz - 1, 1 : Nr - 1] / r[1 : Nr - 1])
                * (T_k[Nz - 1, 2:Nr] - T_k[Nz - 1, 0 : Nr - 2])
                / (2 * dr)
                + (k_eff[Nz - 1, 2:Nr] - k_eff[Nz - 1, 0 : Nr - 2])
                * (T_k[Nz - 1, 2:Nr] - T_k[Nz - 1, 0 : Nr - 2])
                / (4 * dr**2)
                + (k_eff[Nz - 1, 1 : Nr - 1] - k_eff[Nz - 2, 1 : Nr - 1])
                * (T_top[1 : Nr - 1] - T_k[Nz - 2, 1 : Nr - 1])
                / (4 * dz**2)
                + k_eff[Nz - 1, 1 : Nr - 1]
                * (
                    T_k[Nz - 1, 2:Nr]
                    - 2 * T_k[Nz - 1, 1 : Nr - 1]
                    + T_k[Nz - 1, 0 : Nr - 2]
                )
                / dr**2
                + k_eff[Nz - 1, 1 : Nr - 1]
                * (
                    T_top[1 : Nr - 1]
                    - 2 * T_k[Nz - 1, 1 : Nr - 1]
                    + T_k[Nz - 2, 1 : Nr - 1]
                )
                / dz**2
            ) * BETA[
                Nz - 1, 1 : Nr - 1
            ] ** (
                -1
            )

            # Outer boundary points: r = R.
            # curved cylinder surface in contact with air: 0 < z < H
            T_new[1 : Nz - 1, Nr - 1] = T_k[1 : Nz - 1, Nr - 1] + (
                dt / (cp_eff[1 : Nz - 1, Nr - 1] * rho_l)
            ) * (
                (k_eff[1 : Nz - 1, Nr - 1] / r[Nr - 1])
                * (T_edge[1 : Nz - 1] - T_k[1 : Nz - 1, Nr - 2])
                / (2 * dr)
                + (k_eff[1 : Nz - 1, Nr - 1] - k_eff[1 : Nz - 1, Nr - 2])
                * (T_edge[1 : Nz - 1] - T_k[1 : Nz - 1, Nr - 2])
                / (4 * dr**2)
                + (k_eff[2:Nz, Nr - 1] - k_eff[0 : Nz - 2, Nr - 1])
                * (T_k[2:Nz, Nr - 1] - T_k[0 : Nz - 2, Nr - 1])
                / (4 * dz**2)
                + k_eff[1 : Nz - 1, Nr - 1]
                * (
                    T_edge[1 : Nz - 1]
                    - 2 * T_k[1 : Nz - 1, Nr - 1]
                    + T_k[1 : Nz - 1, Nr - 2]
                )
                / dr**2
                + k_eff[1 : Nz - 1, Nr - 1]
                * (
                    T_k[2:Nz, Nr - 1]
                    - 2 * T_k[1 : Nz - 1, Nr - 1]
                    + T_k[0 : Nz - 2, Nr - 1]
                )
                / dz**2
            ) * BETA[
                1 : Nz - 1, Nr - 1
            ] ** (
                -1
            )

            # Center boundary points: r = R.
            # symmetry boundary condition in the centre of the cylinder: 0 < z < H
            T_new[1 : Nz - 1, 0] = T_k[1 : Nz - 1, 0] + (
                dt / (cp_eff[1 : Nz - 1, 0] * rho_l)
            ) * (
                2
                * k_eff[1 : Nz - 1, 0]
                * (T_k[1 : Nz - 1, 1] - 2 * T_k[1 : Nz - 1, 0] + T_center[1 : Nz - 1])
                / dr**2
                + (k_eff[1 : Nz - 1, 1] - k_eff[1 : Nz - 1, 0])
                * (T_k[1 : Nz - 1, 1] - T_center[1 : Nz - 1])
                / (4 * dr**2)
                + (k_eff[2:Nz, 0] - k_eff[0 : Nz - 2, 0])
                * (T_k[2:Nz, 0] - T_k[0 : Nz - 2, 0])
                / (4 * dz**2)
                + k_eff[1 : Nz - 1, 0]
                * (T_k[2:Nz, 0] - 2 * T_k[1 : Nz - 1, 0] + T_k[0 : Nz - 2, 0])
                / dz**2
            ) * BETA[
                1 : Nz - 1, 0
            ] ** (
                -1
            )

            # Bulk of the vial points: : 0 < z < H & 0 < r < R.
            # all the other points in the domain
            T_new[1 : Nz - 1, 1 : Nr - 1] = T_k[1 : Nz - 1, 1 : Nr - 1] + (
                dt / (cp_eff[1 : Nz - 1, 1 : Nr - 1] * rho_l)
            ) * (
                (k_eff[1 : Nz - 1, 1 : Nr - 1] / r[1 : Nr - 1])
                * (T_k[1 : Nz - 1, 2:Nr] - T_k[1 : Nz - 1, 0 : Nr - 2])
                / (2 * dr)
                + (k_eff[1 : Nz - 1, 2:Nr] - k_eff[1 : Nz - 1, 0 : Nr - 2])
                * (T_k[1 : Nz - 1, 2:Nr] - T_k[1 : Nz - 1, 0 : Nr - 2])
                / (4 * dr**2)
                + (k_eff[2:Nz, 1 : Nr - 1] - k_eff[0 : Nz - 2, 1 : Nr - 1])
                * (T_k[2:Nz, 1 : Nr - 1] - T_k[0 : Nz - 2, 1 : Nr - 1])
                / (4 * dz**2)
                + k_eff[1 : Nz - 1, 1 : Nr - 1]
                * (
                    T_k[1 : Nz - 1, 2:Nr]
                    - 2 * T_k[1 : Nz - 1, 1 : Nr - 1]
                    + T_k[1 : Nz - 1, 0 : Nr - 2]
                )
                / dr**2
                + k_eff[1 : Nz - 1, 1 : Nr - 1]
                * (
                    T_k[2:Nz, 1 : Nr - 1]
                    - 2 * T_k[1 : Nz - 1, 1 : Nr - 1]
                    + T_k[0 : Nz - 2, 1 : Nr - 1]
                )
                / dz**2
            ) * BETA[
                1 : Nz - 1, 1 : Nr - 1
            ] ** (
                -1
            )

            # update temperatures
            T_k = T_new

            # Check for local supercooling.
            LCS_i = T_k < T_eq_l
            # the not so supercool part
            LCS_i_r = ~LCS_i

            # ice mass distribution
            m_ice = (
                np.zeros((Nz, Nr)) * LCS_i_r
                + (mass_water - mass_solute * (k_f / M_s) / (T_m - T_new)) * LCS_i
            )
            # ice mass fraction distribution
            w_i_new = m_ice / (mass_water + mass_solute)

            # sparse saving
            if np.mod(i, np.ceil(Nt_solid / N_save_solid)) == 0:
                # saving temperature profiles
                solid_temp_results[i_save, :, :] = T_new - 273.15
                # saving ice mass fractions
                solid_ice_results[i_save, :, :] = w_i_new
                # shelf temperature
                solid_shelf_results[i_save] = T_shelf - 273.15
                # saving sparse time array
                solid_time_results[i_save] = t_nuc + dt * i
                # update saving index
                i_save = i_save + 1

            # solidification sigma
            w_i_r = 2 * np.pi * simps(r * w_i_new, r)
            w_i = simps(w_i_r, z)
            sigma_new = (
                (1 / (np.pi * r[-1] ** 2 * z[-1])) * w_i * mass / (mass - mass_solute)
            )

            # save solidification time
            if (self._stats["t_sol"] is None) & (sigma_new >= 0.9):
                Nt_sol_end = i
                self._stats["t_sol"] = dt * i / 60
                self._stats["t_fr"] = (t_nuc + dt * i) / 60

        # check if solidification is completed
        if Nt_sol_end is None:
            raise ValueError("Solidification is not completed, prolong process time.")

        # temperature distributions in the product
        temp_results = np.concatenate(
            (
                cooling_temp_results[: i_save_end + 1, :],
                solid_temp_results[: i_save - 1, :],
            ),
            axis=0,
        )

        # ice mass fraction results
        ice_results = np.concatenate(
            (
                cooling_ice_results[: i_save_end + 1, :],
                solid_ice_results[: i_save - 1, :],
            ),
            axis=0,
        )

        # total time array
        time_results = (
            np.concatenate(
                (
                    cooling_time_results[: i_save_end + 1],
                    solid_time_results[: i_save - 1],
                ),
                axis=0,
            )
            / 3600
        )

        # total shelf temperature array
        shelf_results = np.concatenate(
            (
                cooling_shelf_results[: i_save_end + 1],
                solid_shelf_results[: i_save - 1],
            ),
            axis=0,
        )

        # save reuslts
        self._domain = (z, r)
        self._prop = (T_eq_l, w_s)
        self._nuc = (t_nuc, T_nuc)
        self._time = time_results
        self._shelfTemp = shelf_results
        self._temp = temp_results
        self._iceMassFraction = ice_results

        return [
            self._stats["T_nuc_min"],
            self._stats["T_nuc_kin"],
            self._stats["T_nuc_mean"],
            self._stats["T_nuc_max"],
            self._stats["t_nuc"],
            self._stats["t_sol"],
            self._stats["t_fr"],
        ]

    # -------------------------------------------------------------------------- #
    # Plotting functions
    # -------------------------------------------------------------------------- #

    # plot distribution of the desired variable
    def plot_cdf(self, what: str = "t_nuc", comp: bool = False):
        """Create CDF plots for Snowing object.

        Args:
            what (str, optional): What to plot, i.e., keys of the stats dict.
                Valid options are t_nuc, T_nuc, t_sol, t_fr.
                Defaults to "t_nuc". Also accepting t_nucleation, T_nucleation
                and t_solidification.
            comp (bool, optional): If complementary CDF is plotted. Defaults to False.

        Args:
            what (str, optional): What to plot, i.e., keys of the stats dict.
                Valid options are t_nuc, T_nuc, t_sol, t_fr.
                Defaults to "t_nuc". Also accepting t_nucleation, T_nucleation
                and t_solidification.

        Raises:
            NotImplementedError: This plot not available for single vial simulation.
        """

        # plot settings
        plt.rcParams.update(
            {
                "axes.axisbelow": True,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.bottom": True,
                "xtick.top": True,
                "ytick.left": True,
                "ytick.right": True,
            }
        )

        # enable this plot only for Nrep > 1
        if self.Nrep <= 1:
            raise NotImplementedError(
                "This plot only available for multiple simulations. "
                + "Run with Nrep > 1."
            )

        # plot settings
        plt.rcParams.update(
            {
                "lines.linewidth": 3,
                "axes.axisbelow": True,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.bottom": True,
                "xtick.top": True,
                "ytick.left": True,
                "ytick.right": True,
            }
        )

        # align the names with Snowfall and Snowflake
        if what == "t_nucleation":
            what = "t_nuc"
        elif what == "T_nucleation":
            what = "T_nuc"
        elif what == "t_solidification":
            what = "t_sol"

        # if nucleation temperatures requested change complementary to True
        if what == "T_nuc":
            comp = True

        # get data
        self._statsMultiple_df = self.results

        # plot attributes
        plt.figure()
        if what == "T_nuc":
            plt.xlabel("nucleation temperature, $T^{nuc}$ [°C]")
            plt.ylabel("nucleation temperature CDF, $F_{nT}$ [K$^{-1}$]")
        elif what == "t_nuc":
            plt.xlabel("nucleation time, $t^{nuc}$ [min]")
            plt.ylabel("nucleation time CDF, $F_{nt}$ [min$^{-1}$]")
        elif what == "t_sol":
            plt.xlabel("solidification time, $t^{sol}$ [min]")
            plt.ylabel("solidification time CDF, $F_{sol}$ [min$^{-1}$]")
        elif what == "t_fr":
            plt.xlabel("complete freezing time, $t^{fr}$ [min]")
            plt.ylabel("frrozen product time CDF, $F_{fr}$ [min$^{-1}$]")
        else:
            raise ValueError(
                'Property to plot incorrectly specified. Use "T_nuc", "t_nuc", "t_sol" or "t_fr" instead.'
            )

        # plot different nucleation temperatures for spatial models
        if (self.const["dimensionality"] != "homogeneous") & (what == "T_nuc"):
            self._data_df_Tnuc = self._statsMultiple_df[
                ["T_nuc_min", "T_nuc_kin", "T_nuc_mean", "T_nuc_max"]
            ]
            sns.ecdfplot(data=self._data_df_Tnuc, complementary=comp)
        # in other cases, just plot the desired quantity
        else:
            sns.ecdfplot(data=self._statsMultiple_df, x=what, complementary=comp)
        plt.grid(axis="both", color="lightgray", linestyle="solid")
        plt.ylim(-0.05, 1.05)
        plt.show()

    # plot the stats (categorical plots) of whichever desired variable
    def plot(self, what: str = "t_nuc", kind: str = "box"):
        """Create categorical plots for Snowing object.

        Args:
            what (str, optional): What to plot, i.e., keys of the stats dict.
                Valid options are t_nuc, T_nuc, t_sol, t_fr.
                Defaults to "t_nuc". Also accepting t_nucleation, T_nucleation
                and t_solidification.
            kind (str, optional): Any sns.catplot 'kind' input is allowed.
                Defaults to "box".
        """

        # enable this plot only for Nrep > 1
        if self.Nrep <= 1:
            raise NotImplementedError(
                "This plot is only available for multiple simulations. "
                + "Run with Nrep > 1."
            )

        # plot settings
        plt.rcParams.update(
            {
                "axes.axisbelow": True,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.bottom": True,
                "xtick.top": True,
                "ytick.left": True,
                "ytick.right": True,
            }
        )

        # align the names with Snowfall and Snowflake
        if what == "t_nucleation":
            what = "t_nuc"
        elif what == "T_nucleation":
            what = "T_nuc"
        elif what == "t_solidification":
            what = "t_sol"

        # get data
        self._statsMultiple_df = self.results

        # plot different nucleation temperatures for spatial models
        if (self.const["dimensionality"] != "homogeneous") & (what == "T_nuc"):
            self._data_df_Tnuc = self._statsMultiple_df[
                ["T_nuc_min", "T_nuc_kin", "T_nuc_mean", "T_nuc_max"]
            ]
            sns.catplot(data=self._data_df_Tnuc, kind=kind)
        # in other cases, just plot the desired quantity
        else:
            sns.catplot(data=self._statsMultiple_df, x=what, kind=kind)

    # plot the temperature evolution
    def _plot_temperature_evolution(self):
        # temperature evolution plot
        plt.figure()
        # plot nucleation time and equilibrum temperature
        plt.plot(
            [self._stats["t_nuc"] / 60, self._stats["t_nuc"] / 60],
            [self._shelfTemp.min() - 5, self._shelfTemp.max() + 5],
            "--",
            color="red",
            linewidth=1,
            label="$t = t^{nuc}$",
        )
        # plot time of frozen product
        if self._stats["t_fr"] is not None:
            plt.plot(
                [self._stats["t_fr"] / 60, self._stats["t_fr"] / 60],
                [self._shelfTemp.min() - 5, self._shelfTemp.max() + 5],
                "--",
                color="indigo",
                linewidth=1,
                label="$t = t^{fr}$",
            )
        plt.plot(
            [-0.05 * self._time.max(), 1.05 * self._time.max()],
            [self._prop[0] - 273.15, self._prop[0] - 273.15],
            "--",
            color="magenta",
            linewidth=1,
            label="$T = T^{eq}_{l}$",
        )
        # plot shelf temperature evolution
        plt.plot(self._time, self._shelfTemp, "b-", linewidth=2.5, label="shelf T")

        if self.const["dimensionality"] == "homogeneous":
            # plot temperature evolution for homogeneous model
            plt.plot(self._time, self._temp, "k-", linewidth=2.5, label="product T")

        # access legend objects automatically created from data
        handles, labels = plt.gca().get_legend_handles_labels()

        # plot product temperature evolution for 1D case
        if self.const["dimensionality"] == "spatial_1D":
            # colors
            colors, my_cmap, norm = Utils.colormap(self._domain)
            # plot temperature evolution for different vertical positions
            for i in range(0, len(self._domain), round(len(self._domain) / 11)):
                plt.plot(self._time, self._temp[:, i], color=colors[i], linewidth=2.5)
            # colorbar
            sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            cbar = plt.colorbar(sm, ax=plt.gca(), aspect=35)
            cbar.set_label("dimensionless vertical position, z/H [-]", rotation=90)
            # add colormap to the legend entries
            cmaps_gradients = my_cmap(np.linspace(0, 1, 30))
            # combine patches into a list and add to the handles and labels
            cmap_patches = [
                patches.Patch(facecolor=c, edgecolor=c) for c in cmaps_gradients
            ]
            handles.append(cmap_patches)
            labels.append("product T")

        # plot product temperature evolution for 2D case
        elif self.const["dimensionality"] == "spatial_2D":
            # colors
            colors, my_cmap, norm = Utils.colormap(self._domain[0])
            # plot temperature evolution for different vertical positions
            for i in range(0, len(self._domain[0]), round(len(self._domain[0]) / 11)):
                plt.plot(
                    self._time, self._temp[:, i, -1], color=colors[i], linewidth=2.5
                )
            # colorbar
            sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            cbar = plt.colorbar(sm, ax=plt.gca(), aspect=35)
            cbar.set_label(
                "dimensionless vertical position (at r = 0), z/H [-]", rotation=90
            )
            # add colormap to the legend entries
            cmaps_gradients = my_cmap(np.linspace(0, 1, 30))
            # combine patches into a list and add to the handles and labels
            cmap_patches = [
                patches.Patch(facecolor=c, edgecolor=c) for c in cmaps_gradients
            ]
            handles.append(cmap_patches)
            labels.append("product T")

        # plot labels and attributes
        plt.xlabel("process time, $t$ [h]")
        plt.ylabel("temperature, $T$ [°C]", labelpad=8)
        plt.grid(axis="both", color="lightgray", linestyle="solid")
        plt.xlim(-0.05 * self._time.max(), 1.05 * self._time.max())
        plt.ylim(self._shelfTemp.min() - 5, self._shelfTemp.max() + 5)
        plt.title("Product temperature evolution")
        plt.legend(
            handles=handles,
            labels=labels,
            handler_map={list: HandlerTuple(ndivide=None, pad=0)},
            loc="best",
        )

    # plot ice mass fraction evolution
    def _plot_ice_mass_fraction_evolution(self):
        # ice mass fraction plot
        plt.figure()
        # plot lines for nucleation time and maximum ice mass fraction
        plt.plot(
            [self._stats["t_nuc"] / 60, self._stats["t_nuc"] / 60],
            [-0.05, 1.05],
            "--",
            color="red",
            linewidth=1,
            label="$t = t^{nuc}$",
        )
        # plot time of frozen product
        if self._stats["t_fr"] is not None:
            plt.plot(
                [self._stats["t_fr"] / 60, self._stats["t_fr"] / 60],
                [-0.05, 1.05],
                "--",
                color="indigo",
                linewidth=1,
                label="$t = t^{fr}$",
            )
        plt.plot(
            [-0.05 * self._time.max(), 1.05 * self._time.max()],
            [1 - self._prop[1], 1 - self._prop[1]],
            "--",
            color="blue",
            linewidth=1,
            label="$w_{i} = 1 - w_{s}$",
        )

        if self.const["dimensionality"] == "homogeneous":
            # plot temperature evolution for homogeneous model
            plt.plot(
                self._time,
                self._iceMassFraction,
                "k-",
                linewidth=2.5,
                label="product $w_{i}$",
            )

        # access legend objects automatically created from data
        handles, labels = plt.gca().get_legend_handles_labels()

        # plot product temperature evolution for 1D case
        if self.const["dimensionality"] == "spatial_1D":
            # colors
            colors, my_cmap, norm = Utils.colormap(self._domain)
            # colorbar
            sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            cbar = plt.colorbar(sm, ax=plt.gca(), aspect=35)
            cbar.set_label("dimensionless vertical position, z/H [-]", rotation=90)
            # plot temperature evolution for different vertical positions
            for i in range(0, len(self._domain), round(len(self._domain) / 11)):
                plt.plot(
                    self._time,
                    self._iceMassFraction[:, i],
                    color=colors[i],
                    linewidth=2.5,
                )
            # add colormap to the legend entries
            cmaps_gradients = my_cmap(np.linspace(0, 1, 30))
            # combine patches into a list and add to the handles and labels
            cmap_patches = [
                patches.Patch(facecolor=c, edgecolor=c) for c in cmaps_gradients
            ]
            handles.append(cmap_patches)
            labels.append("product $w_{i}$")

        # plot product temperature evolution for 2D case
        elif self.const["dimensionality"] == "spatial_2D":
            # colors
            colors, my_cmap, norm = Utils.colormap(self._domain[0])
            # colorbar
            sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            cbar = plt.colorbar(sm, ax=plt.gca(), aspect=35)
            cbar.set_label(
                "dimensionless vertical position (at r = 0), z/H [-]", rotation=90
            )
            # plot temperature evolution for different vertical positions
            for i in range(0, len(self._domain[0]), round(len(self._domain[0]) / 11)):
                plt.plot(
                    self._time,
                    self._iceMassFraction[:, i, -1],
                    color=colors[i],
                    linewidth=2.5,
                )
            # add colormap to the legend entries
            cmaps_gradients = my_cmap(np.linspace(0, 1, 30))
            # combine patches into a list and add to the handles and labels
            cmap_patches = [
                patches.Patch(facecolor=c, edgecolor=c) for c in cmaps_gradients
            ]
            handles.append(cmap_patches)
            labels.append("product $w_{i}$")

        # plot labels and attributes
        plt.xlabel("process time, $t$ [h]")
        plt.ylabel("ice mass fraction, $w_i$ [-]")
        plt.grid(axis="both", color="lightgray", linestyle="solid")
        plt.xlim(-0.05 * self._time.max(), 1.05 * self._time.max())
        plt.ylim(-0.05, 1.05)
        plt.title("Ice mass fraction evolution")
        plt.legend(
            handles=handles,
            labels=labels,
            handler_map={list: HandlerTuple(ndivide=None, pad=0)},
            loc="best",
        )

    # plot evolution function
    def plot_evolution(self, what: str = "temperature"):
        """Function to plot the spatial evolution of temperature or ice mass fractions.

        Args:
            what (str, optional): Quantity to be plotted. Defaults to "temperature".

        Raises:
            NotImplementedError: This plotting functionality is only implemented for single simulations (when Nrep = 1), otherwise, an error is raised.
            ValueError: Raises error if quantity for plotting is incorrectly specified.
        """

        # plot settings
        plt.rcParams.update(
            {
                "axes.axisbelow": True,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.bottom": True,
                "xtick.top": True,
                "ytick.left": True,
                "ytick.right": True,
            }
        )

        if self.Nrep > 1:
            raise NotImplementedError(
                "This plot only available for single simulations. "
                + "Run with Nrep = 1."
            )

        if what == "temperature":
            self._plot_temperature_evolution()
        elif what == "ice_mass_fraction":
            self._plot_ice_mass_fraction_evolution()
        else:
            raise ValueError(
                'Property to plot incorrectly specified. Use "temperature" or "ice_mass_fraction" instead.'
            )

    # repr function
    def __repr__(self) -> str:
        """Return string representation of the Snowing class.

        Returns:
            str: The Snowing class string representation giving some basic info.
        """
        return (
            f"Snowing([{self.Nrep} Snowing{'s' if self.Nrep > 1 else ''}, "
            + f"dimensionality: {self.const['dimensionality']}, "
            + f"configuration: {self.const['configuration']}])"
        )
