"""Implement Snowflake class.

This module contains the Snowflake class used to run simulations
of water nucleation in vials.
"""

from ethz_snow.operatingConditions import OperatingConditions
from ethz_snow.constants import calculateDerived

import numpy as np
import pandas as pd
import re
import warnings

# import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from scipy.sparse import csr_matrix
from typing import List, Tuple, Type, Union, Sequence, Optional

from ethz_snow.__init__ import __citation__, __bibtex__

HEATFLOW_REQUIREDKEYS = ("int", "ext", "s0")
VIAL_GROUPS = (
    "corner",
    "edge",
    "core",
    "side",
    "all",
    "center",
)  # center will be deprecated


class Snowflake:
    """A class to handle a single Stochastic Nucleation of Water simulation.

    More information regarding the equations and their derivation can be found in
    "Stochastic shelf-scale modeling framework for the freezing stage in
    freeze-drying processes", Deck, Ochsenbein, and Mazzotti (2022),
    Int J Pharm, 613, 121276, https://doi.org/10.1016/j.ijpharm.2021.121276.

    Parameters:
        configPath (Optional[str]): The path of the (optional) custom config yaml.
        const (dict): A dictionary of constants to be used.
        dt (float): Time step.
        H_ext (np.ndarray): External heat transfer vector.
        H_int (csr_matrix): Internal heat transfer matrix.
        H_shelf (np.ndarray): Shelf heat transfer vector.
        initIce (str, optional): Which formulation for the initial amount of ice
                formed to use (see docs). Can be 'direct' or 'indirect'.
                Defaults to 'indirect'.
        k (dict): Heat transfer coefficients.
        N_vials (tuple): Number of vials in each dimension.
        N_vials_total (int): Total number of vials.
        opcond (OperatingConditions): Operating conditions of run.
        seed (int): Seed to be used in rng.
        seed_v (int): Seed to be used in rng for vial-dependent parameters.
        simulationStatus (int): Status of simulation (0 = not run, 1 = run).
        solidificationThreshold (float): What sigma value constitutes a 'solid'.
        stats (dict): Run statistics (nucleation time, etc.).
        T_k_0 (float): Initial temperature of vials.
        X_sigma (np.ndarray): Sigma state over time.
        X_T (np.ndarray): Temperature state over time.
    """

    def __init__(
        self,
        k: dict = {"int": 20, "ext": 20, "s0": 20, "s_sigma_rel": 0.1},
        N_vials: Tuple[int, int, int] = (
            7,
            7,
            1,
        ),  # should this be part of operating conditions?
        initialStates: Optional[dict] = {"temp": None, "sigma": None},
        storeStates: Union[str, Sequence[str], Sequence[int], None] = None,
        solidificationThreshold: float = 0.9,
        dt: float = 2,
        seed: int = 2021,
        seed_v: int = 2024,
        opcond: OperatingConditions = OperatingConditions(),
        configPath: Optional[str] = None,
        initIce: str = "indirect",
    ):
        """Construct a Snowflake object.

        Args:
            k (dict, optional): A dictionary containing the heat transfer coefficients.
                Must contain keys 'int', 'ext', 's0'.
                Defaults to {"int": 20, "ext": 20, "s0": 20, "s_sigma_rel": 0.1}.
            N_vials (Tuple[int, int, int], optional): Number of vials in each dimension.
                Defaults to ( 7, 7, 1).
            initialStates (dict): Initial states of the vials (temp and sigma).
            storeStates (Union[str, Sequence[str], Sequence[int], None]):
                Whether and for which vials to store states (temp, sigma).
                Defaults to None -> Nothing stored.
            solidificationThreshold (float, optional): What sigma value is
                a 'solid'. Defaults to 0.9.
            dt (float, optional): Time step size. Defaults to 2.
            seed (int, optional): Seed to use for rng. Defaults to 2021.
            seed_v (int, optional): Vial-dependent seed to use for rng. Defaults to 2024.
            opcond (OperatingConditions, optional): Operating conditions to apply.
                Defaults to OperatingConditions().
            configPath (Optional[str], optional): The location of a custom configuration
                file (in yaml format), used to update/overwrite the default settings.
            initIce (str, optional): Which formulation for the initial amount of ice
                formed to use (see docs). Defaults to 'indirect'.

        Raises:
            TypeError: If k is not a dict.
            ValueError: If k does not contain all necessary keys.
            ValueERror: If len(N_vials) is not 3.
            TypeError: If opcond is not of type operatingConditions.
            ValueError: If initialstates is malformed.
            NotImplementedError: If initialstates has non-scalar temp or non-None sigma.
            ValueError: If storeStates is not meaningful.

        Examples:
            >>> Sf = Snowflake()
            >>> Sf = Snowflake(Nvials=(4, 3, 1), dt = 5)
            >>> Sf = Snowflake(storeStates = ['edge_random_4', 'uniform.core.5'])
            >>> Sf = Snowflake(storeStates = (0, 4, 10))
        """

        if isinstance(N_vials, list):
            # we exepct this to be immutable
            N_vials = tuple(N_vials)

        if len(N_vials) != 3:
            raise ValueError(f"Length of N_vials must be 3, was {len(N_vials)}.")
        else:
            self.N_vials = N_vials

        if self.N_vials[2] > 1:
            # 3D case, k_shelf can be missing
            hflowKeys = tuple([elem for elem in HEATFLOW_REQUIREDKEYS if elem != "s0"])
        else:
            hflowKeys = HEATFLOW_REQUIREDKEYS

        if not isinstance(k, dict):
            raise TypeError(f"Input k must be of type dict. Was given {type(k)}.")
        elif not all([key in k.keys() for key in hflowKeys]):
            raise ValueError(
                (
                    f"A required key was missing from dictionary k, specifically "
                    + f"{set(hflowKeys) - set(k.keys())}."
                )
            )
        else:
            self.k = k

        self.dt = dt

        if not isinstance(opcond, OperatingConditions):
            raise TypeError(
                "Input opcond must be an instance of class OperatingConditions."
            )
        else:
            self.opcond = opcond

        self.solidificationThreshold = solidificationThreshold
        self._X = None
        self._t = None
        self.stats = dict()

        self._simulationStatus = 0

        self._H_int = None
        self._H_ext = None
        self._H_shelf = None

        self.configPath = configPath

        # store the seed to look it up if need be
        self.seed = seed
        self.seed_v = seed_v

        # remember what N_vials was used to build heat flow matrices
        # so if it changes we know to rebuild them
        self._NvialsUsed = self.N_vials

        if (initialStates is not None) and (
            (not isinstance(initialStates, dict))
            or ("temp" not in initialStates.keys())
        ):
            raise ValueError(
                "initialStates malformed. Must be None or dict with keys "
                + f"'temp' (and 'sigma'). Was {initialStates}."
            )

        if (initialStates is None) or (initialStates["temp"] is None):
            self.T_k_0 = self.opcond.cooling["start"]  # C
        elif not hasattr(initialStates["temp"], "__len__"):
            self.T_k_0 = initialStates["temp"]  # C
        elif (hasattr(initialStates["temp"], "__len__")) and (
            len(initialStates["temp"]) == 1
        ):
            self.T_k_0 = initialStates["temp"][0]
        elif (hasattr(initialStates["temp"], "__len__")) and (
            len(initialStates["temp"]) > 1
        ):
            raise NotImplementedError(
                "Currently can only deal with constant initial temp."
            )

        if ("sigma" in initialStates.keys()) and (initialStates["sigma"] is not None):
            raise NotImplementedError("Non-zero sigma initial state not implemented.")

        self._emptyStore = False
        if isinstance(storeStates, (list, tuple)):
            if all([isinstance(x, int) for x in storeStates]):
                if any(np.array(storeStates) > self.N_vials_total - 1) or any(
                    np.array(storeStates) < 0
                ):
                    raise ValueError(
                        f"Entries in storeStates must be >0 "
                        + f"and <{self.N_vials_total - 1}."
                    )
                storageMask = np.zeros(self.N_vials_total, dtype=bool)
                storageMask[list(storeStates)] = True
            elif all([isinstance(x, str) for x in storeStates]):
                storageMasks = list(
                    map(lambda x: self._interpretStorageString(x), storeStates)
                )
                storageMask = np.logical_or.reduce(storageMasks)
            else:
                raise ValueError(
                    "storeStates must be a sequence of all int or all str."
                )
        elif isinstance(storeStates, str):
            storageMask = self._interpretStorageString(storeStates)
        elif storeStates is None:
            storageMask = np.zeros(self.N_vials_total, dtype=bool)
            self._emptyStore = True

        self._storageMask = storageMask

        if not isinstance(initIce, str) or (
            (initIce.lower() != "direct") and (initIce.lower() != "indirect")
        ):
            raise ValueError(
                "initIce must be a string equal to 'indirect' or 'direct'. ",
                f"But it was {initIce}.",
            )
        self.initIce = initIce.lower()

    @property
    def simulationStatus(self) -> int:
        """Return simulation status of instance.

        Returns:
            int: Simulation status. 0 = not run, 1 = run.
        """
        if (self._simulationStatus == 1) or (self._X is not None):
            self._simulationStatus = 1

        return self._simulationStatus

    @property
    def H_int(self) -> csr_matrix:
        """Return the internal heat transfer matrix.

        Returns:
            csr_matrix: The internal heat transfer matrix.
        """
        if (self._H_int is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
        return self._H_int

    @property
    def H_ext(self) -> np.ndarray:
        """Return the external heat transfer vector.

        Returns:
            np.ndarray: The external heat transfer vector.
        """
        if (self._H_ext is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
        return self._H_ext

    @property
    def H_shelf(self) -> np.ndarray:
        """Return the shelf heat transfer vector.

        Returns:
            np.ndarray: The shelf heat transfer vector.
        """
        if (
            (self._H_shelf is None)
            or (self.N_vials != self._NvialsUsed)
            or (self.seed != self._seedUsed)
        ):
            self._buildShelfHeatFlow()
        return self._H_shelf

    @property
    def N_vials_total(self) -> int:
        """Return total number of vials in system."""
        return int(np.prod(self.N_vials))

    @property
    def X_T(self) -> np.ndarray:
        """Get the temperature states.

        Returns:
            np.ndarray: An array containing the vial temperatures over time.
                This is just a slice of _X!
        """
        return self._X[: int(np.sum(self._storageMask)), :]

    @property
    def X_sigma(self) -> np.ndarray:
        """Get the sigma states.

        Returns:
            np.ndarray: An array containing the vial sigmas over time.
                This is just a slice of _X!
        """
        return self._X[int(np.sum(self._storageMask)) :, :]

    @property
    def seed(self) -> int:
        """Get or set the random seed.

        Setting the seed value will initialize a new rng under the hood.

        Returns:
            int: The seed of the Snowflake.
        """
        return self._seed

    @seed.setter
    def seed(self, value):
        # best practice is now to store the rng in an object and pass it around
        # https://towardsdatascience.com/stop-using-numpy-random-seed-581a9972805f
        self._rng = np.random.default_rng(value)
        # store the seed to look it up if need be
        self._seed = value

        # let H_shelf figure out whether it needs to be updated
        _ = self.H_shelf

    @property
    def configPath(self) -> Optional[str]:
        """Get configPath property."""
        return self._configPath

    @configPath.setter
    def configPath(self, value):
        """Set configPath property."""
        self.const = calculateDerived(value)
        self._configPath = value

    def _interpretStorageString(self, myString):
        """Interpret the storeState string.

        Returns the interpretation of the storeState string
        indicating which states should be stored in _X.
        """
        myString = myString.lower()
        if not any(
            [word in myString for word in list(VIAL_GROUPS) + ["random", "uniform"]]
        ):
            raise ValueError("No valid vial group or key word used in storeStates.")
        elif any([word in myString for word in list(VIAL_GROUPS)]):
            for group in VIAL_GROUPS:
                if group in myString:
                    storageMask = self.getVialGroup(group)
                    break
        elif ("random" in myString) or ("uniform" in myString):
            # assume vial group 'all' is implied
            storageMask = np.ones(self.N_vials_total, dtype=bool)

        if ("random" in myString) or ("uniform" in myString):
            howMany = re.findall(r"\d+", myString)
            if len(howMany) == 0:
                # unclear how many to pick, default is 10%
                howMany = int(np.ceil(0.1 * self.N_vials_total))
            elif len(howMany) == 1:
                howMany = int(howMany[0])
            elif len(howMany) > 1:
                raise ValueError("storeStates strings must contain one number at most.")

            I_candidates = np.where(storageMask)[0]
            if "random" in myString:
                I_toStore = self._rng.choice(I_candidates, size=howMany, replace=False)
            elif "uniform" in myString:
                stepSize = int(np.ceil(len(I_candidates) / howMany))
                I_fromCandidates = np.arange(0, len(I_candidates), stepSize, dtype=int)
                I_toStore = I_candidates[I_fromCandidates]
            storageMask = np.zeros(self.N_vials_total, dtype=bool)
            storageMask[I_toStore] = True

        return storageMask

    def run(self):
        """Run the simulation."""
        # clean up any potential old simulations
        self._X = None
        self._t = None

        # Obtain external temperature profile over time
        T_shelf = self.opcond.tempProfile(self.dt)
        T_ext = T_shelf  # need to make a switch so that this can be decoupled - DRO XX

        # store stuff in local variables to reduce getter method calls
        # note that getter method handles the construction of the matrices
        # in case they haven't been initialized yet
        H_int = self.H_int
        H_ext = self.H_ext
        H_shelf = self.H_shelf

        # constants
        solid_fraction = self.const["solid_fraction"]
        cp_s = self.const["cp_s"]
        cp_w = self.const["cp_w"]
        cp_i = self.const["cp_i"]
        cp_solution = self.const["cp_solution"]
        depression = self.const["depression"]
        mass = self.const["mass"]
        alpha = self.const["alpha"]
        beta_solution = self.const["beta_solution"]
        T_eq = self.const["T_eq"]
        T_eq_l = self.const["T_eq_l"]
        # T_init = self.const["T_init"]
        hl = self.const["hl"]
        # kb = self.const["kb"]
        a = self.const["a"]
        b = self.const["b"]
        c = self.const["c"]
        V = self.const["V"]

        # random numbers to determine vial-dependent nucleation parameters
        np.random.seed(self.seed_v)
        vial_random_numbers = np.random.rand(self.N_vials_total)
        xi_v = norm.ppf(vial_random_numbers)

        # pre-exponential nucleation parameter
        kb = 10 ** (-(a + xi_v * c))

        N_timeSteps = (
            int(np.ceil(self.opcond.t_tot / self.dt)) + 1
        )  # the total number of timesteps
        N_vials_total = self.N_vials_total

        # Initial state of the system
        # initial temperature
        T_k = np.ones(N_vials_total) * self.T_k_0  # [C]
        sigma_k = np.zeros(N_vials_total)  # (0 = liq, 1 = completely frozen)

        # Pre-allocate memory
        # our state matrix containing the entire history of the system
        X = np.zeros((2 * np.sum(self._storageMask), N_timeSteps))
        t = np.arange(0, N_timeSteps * self.dt, self.dt)
        stats = dict(
            t_nucleation=np.full(N_vials_total, np.nan),
            T_nucleation=np.full(N_vials_total, np.nan),
            t_solidification=np.full(N_vials_total, np.nan),
        )

        if np.any(t >= self.opcond.cnt):
            k_CN = np.argmax(t >= self.opcond.cnt)
        else:
            k_CN = N_timeSteps + 1

        # Iterate over time steps
        for k in np.arange(N_timeSteps):

            T_shelf_k = T_shelf[k]
            T_ext_k = T_ext[k]
            # toc1_l = time.perf_counter()
            if not self._emptyStore:
                # does not give a huge boost.
                x_k = np.concatenate([T_k, sigma_k])
                X[:, k] = x_k[np.concatenate([self._storageMask, self._storageMask])]

            # calculate heatflows for all vials
            q_k = (
                H_int @ T_k + H_ext * (T_ext_k - T_k) + H_shelf * (T_shelf_k - T_k)
            )  # [XXX] units? - DRO

            # a mask that is True where a vial is still fully liquid
            liquidMask = sigma_k == 0
            solidMask = ~liquidMask  # for convenience
            # SOLID(IFYING) VIALS
            if any(solidMask):
                solidifiedMask = (sigma_k > self.solidificationThreshold) & np.isnan(
                    stats["t_solidification"]
                )
                stats["t_solidification"][solidifiedMask] = (
                    t[k] - stats["t_nucleation"][solidifiedMask]
                )
                # heat capacity of solidifying vials
                # a vector of cps for all solidifying vials
                cp_sigma = solid_fraction * cp_s + (1 - solid_fraction) * (
                    sigma_k[solidMask] * (cp_i - cp_w) + cp_w
                )
                beta = depression * mass * cp_sigma

                deltaSigma = (
                    q_k[solidMask]
                    / (alpha - beta / (1 - sigma_k[solidMask]) ** 2)
                    * self.dt
                )
                sigma_k[solidMask] = sigma_k[solidMask] + deltaSigma
                # assumption is that during solidification T = Teq
                T_k[solidMask] = T_eq - depression * (1 / (1 - sigma_k[solidMask]))

            # LIQUID VIALS
            if any(liquidMask):
                T_k[liquidMask] = T_k[liquidMask] + q_k[liquidMask] / hl * self.dt
                # a mask that is True where the temperature is
                # below the equilibrium temperature
                superCooledMask = T_k < T_eq_l
                # a vial that is both liquid and supercooled can nucleate
                nucleationCandidatesMask = liquidMask & superCooledMask
                # the total number of nucleation candidates
                n_nucleationCandidates = np.sum(nucleationCandidatesMask)

                # nucleation probabilities for candidates
                P = np.zeros(N_vials_total)
                diceRolls = np.zeros(N_vials_total)

                P[nucleationCandidatesMask] = (
                    kb[nucleationCandidatesMask]
                    * V
                    * (T_eq_l - T_k[nucleationCandidatesMask]) ** b
                    * self.dt
                )
                # when we reach the timepoint of controlled nucleation
                # all vials (that thermodynamically can) nucleate
                if k == k_CN:
                    P.fill(1)
                    frac_justCant = 1 - np.sum(nucleationCandidatesMask) / N_vials_total
                    if frac_justCant >= 0.1:
                        print(
                            (
                                f"WARNING: A significant number of vials "
                                + f"({frac_justCant * 100:4.1f}%) "
                                + "were not supercooled when controlled "
                                + "nucleation triggered."
                            )
                        )
                diceRolls[nucleationCandidatesMask] = self._rng.random(
                    n_nucleationCandidates
                )

                # Nucleation
                nucleatedVialsMask = nucleationCandidatesMask & (diceRolls < P)

                stats["t_nucleation"][nucleatedVialsMask] = t[k] + self.dt
                stats["T_nucleation"][nucleatedVialsMask] = T_k[nucleatedVialsMask]

                q0 = (T_eq_l - T_k[nucleatedVialsMask]) * cp_solution * mass

                if self.initIce == "indirect":
                    sigma_k[nucleatedVialsMask] = -q0 / (alpha - beta_solution)
                else:
                    # LTD: referred to as "direct" approach
                    # this is the one used in the derivation of the manuscripts
                    gamma_direct = -alpha / mass / cp_solution
                    B_direct = T_eq - T_k[nucleatedVialsMask] + gamma_direct
                    C_direct = T_k[nucleatedVialsMask] - T_eq + depression

                    sigma_k[nucleatedVialsMask] = (
                        -B_direct + np.sqrt(B_direct**2 + 4 * gamma_direct * C_direct)
                    ) / (-2 * gamma_direct)

                # assumption is that during solidification T = Teq
                T_k[nucleatedVialsMask] = T_eq - depression * (
                    1 / (1 - sigma_k[nucleatedVialsMask])
                )

        self.stats = stats
        self._X = X  # store the state matrix
        self._t = t  # store the time vector

    def nucleationTimes(
        self, group: Union[str, Sequence[str]] = "all", fromStates: bool = False
    ) -> np.ndarray:
        """Return array of nucleation times.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            fromStates (bool, optional): Whether or not to calculate from
                states directly. Defaults to False.

        Returns:
            np.ndarray: The nucleation times.
        """
        if fromStates:
            I_nucleation, lateBloomers = self._sigmaCrossingIndices(threshold=0)

            # nucleation times for all nucleated vials
            # need to make sure this is float so no problems arise later
            t_nucleation_states = self._t[I_nucleation].astype(float)

            # non-nucleated vials are set to NaN
            t_nucleation_states[lateBloomers] = np.nan

            t_nucleation = np.full(self.N_vials_total, np.nan)
            t_nucleation[self._storageMask] = t_nucleation_states
        else:
            t_nucleation = self.stats["t_nucleation"]

        I_groups = self.getVialGroup(group)

        return t_nucleation[I_groups]

    def nucleationTemperatures(
        self, group: Union[str, Sequence[str]] = "all", fromStates: bool = False
    ) -> np.ndarray:
        """Return array of nucleation temperatures.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            fromStates (bool, optional): Whether or not to calculate from
                states directly. Defaults to False.

        Returns:
            np.ndarray: The nucleation temperatures.
        """
        if fromStates:
            I_nucleation, lateBloomers = self._sigmaCrossingIndices(threshold=0)
            T_nucleation_states = np.array(
                [
                    self.X_T[i, I - 1]
                    for i, I in zip(range(self.N_vials_total), I_nucleation)
                ]
            ).astype(
                float
            )  # should always be float, but just to be sure

            # non-nucleated vials are set to NaN
            T_nucleation_states[lateBloomers] = np.nan

            T_nucleation = np.full(self.N_vials_total, np.nan)
            T_nucleation[self._storageMask] = T_nucleation_states

        else:
            T_nucleation = self.stats["T_nucleation"]

        I_groups = self.getVialGroup(group)

        return T_nucleation[I_groups]

    def solidificationTimes(
        self,
        group: Union[str, Sequence[str]] = "all",
        threshold: Optional[float] = None,
        fromStates: bool = False,
    ) -> np.ndarray:
        """Return array of solidification times.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            threshold (Optional[float], optional): The threshold used to
                define 'solidified'. Defaults to self.solidificationThreshold.
            fromStates (bool, optional): Whether or not to calculate from
                states directly. Defaults to False.

        Returns:
            np.ndarray: The solidification times.
        """
        if threshold is None:
            threshold = self.solidificationThreshold

        if fromStates:
            t_nucleation = self.nucleationTimes(fromStates=fromStates)

            I_solidification, neverGrownUp = self._sigmaCrossingIndices(
                threshold=threshold
            )

            # nucleation times for all nucleated vials
            # need to make sure this is float so no problems arise later
            t_solidified = self._t[I_solidification].astype(float)

            # never fully solidified vials are set to NaN
            t_solidified[neverGrownUp] = np.nan

            # solidification is the difference between solidified and nucleated times
            t_solidification_states = t_solidified

            t_solidification = np.full(self.N_vials_total, np.nan)
            t_solidification[self._storageMask] = t_solidification_states

            t_solidification = t_solidification - t_nucleation
        else:
            if threshold != self.solidificationThreshold:
                raise ValueError(
                    "Threshold cannot differ from solidificationThreshold when "
                    + "not calculating from states."
                )
            t_solidification = self.stats["t_solidification"]

        I_groups = self.getVialGroup(group)

        return t_solidification[I_groups]

    def sigmaCounter(
        self,
        time: Union[Sequence[float], float],
        threshold: Optional[float] = None,
        fromStates: bool = False,
    ) -> np.ndarray:
        """Return counter of vials satisfying sigma>threshold at time.

        Args:
            time (Union[Sequence[float], float]): The time(s) for which to
                return the counter.
            threshold (Optional[float], optional): The threshold to apply.
                Defaults to self.solidificationThreshold.
            fromStates (bool, optional): Whether or not to calculate from
                states directly. Defaults to False.

        Raises:
            ValueError: Simulation needs to be run first.
            ValueError: fromStates is False and threshold is
                neither 0 or self.solidificationThreshold.

        Returns:
            np.ndarray: [description]
        """
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run first.")

        if threshold is None:
            # default is to give solidified vials
            threshold = self.solidificationThreshold

        if isinstance(time, (float, int)):
            # we iterate over this,
            # so it must be a list
            time = [time]

        counter = np.zeros(len(time))
        for i, t in enumerate(time):

            if fromStates:
                I_time = np.argmax(self._t >= t)
                counter[i] = np.sum(self.X_sigma[:, I_time] > threshold, axis=0)
            else:
                if (threshold > 0) and (threshold != self.solidificationThreshold):
                    raise ValueError(
                        "Threshold cannot differ from solidificationThreshold when "
                        + "not calculating from states."
                    )
                elif threshold > 0:
                    counter[i] = np.sum(self.stats["t_solidification"] <= t)
                elif threshold == 0:
                    counter[i] = np.sum(self.stats["t_nucleation"] <= t)

        return counter

    def getVialGroup(
        self, group: Union[str, Sequence[str]] = "all"
    ) -> np.ndarray:  # @ DRO: Need to adjust definitions for corner, edge, side etc...
        """Return mask for given group.

        mask[i] = True iff vial[i] is in G where G can be a list of groups.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Can be "corner", "edge", "side", "core", "all". Defaults to "all".

        Raises:
            ValueError: Group is not known.

        Returns:
            np.ndarray: A boolean mask of size (N_vials_tot,).

        Examples:
            >>> S = Snowflake(N_vials = (2, 2, 1))
            >>> S.getVialGroup('corner')
            array([True, True, True, True])

            >>> S.getVialGroup(['edge', 'core'])
            array([False, False, False, False])
        """
        if isinstance(group, str):
            group = [group]

        _, VIAL_EXT = self._buildInteractionMatrices()

        myMask = np.zeros(self.N_vials_total, dtype=bool)

        # square vial arrangement
        if self.const["vial_arrangement"] == "square":
            for g in group:
                if g == "corner":
                    # for the 3D case, a corner vial has 3 external IAs
                    myMask = myMask | (VIAL_EXT == 3 - (self.N_vials[2] == 1))
                elif g == "edge":
                    myMask = myMask | (VIAL_EXT == 2 - (self.N_vials[2] == 1))
                elif g == "side":
                    # side and edge are equivalent for 2D
                    myMask = myMask | (VIAL_EXT == 1)
                elif g == "core":
                    myMask = myMask | (VIAL_EXT == 0)
                elif g == "center":
                    myMask = myMask | (VIAL_EXT == 0)
                    warnings.warn(
                        "'center' is deprecated terminology. Please use 'core' instead.",
                        DeprecationWarning,
                    )
                elif g == "all":
                    myMask = np.ones(self.N_vials_total, dtype=bool)
                    break
                else:
                    raise ValueError(f"Group {g} is not a known vial group.")

        # hexagonal vial arrangement: only considering corner, edge and core
        else:
            for g in group:
                if g == "corner":
                    # corner vials have 4 ext. IAs in 2D batch
                    myMask = myMask | (VIAL_EXT == 5 - (self.N_vials[2] == 1))
                elif g == "edge":
                    myMask = myMask | (
                        (VIAL_EXT > 0) & (VIAL_EXT < 5 - (self.N_vials[2] == 1))
                    )
                elif g == "side":
                    # side and edge are equivalent for the hexagonal
                    myMask = myMask | (
                        (VIAL_EXT > 0) & (VIAL_EXT < 5 - (self.N_vials[2] == 1))
                    )
                elif g == "core":
                    myMask = myMask | (VIAL_EXT == 0)
                elif g == "center":
                    myMask = myMask | (VIAL_EXT == 0)
                    warnings.warn(
                        "'center' is deprecated terminology. Please use 'core' instead.",
                        DeprecationWarning,
                    )
                elif g == "all":
                    myMask = np.ones(self.N_vials_total, dtype=bool)
                    break
                else:
                    raise ValueError(f"Group {g} is not a known vial group.")

        return myMask

    def plot(
        self,
        kind: str = "box",
        what: str = "t_nucleation",
        group: Union[str, Sequence[str]] = "all",
        context: bool = True,
    ):
        """Create plots for Snowflake object.

        Args:
            kind (str, optional): Type of plot to print. If input is
                'trajectories' will print temperature or sigma profile.
                Otherwise will print from stats dict. In that case,
                any sns.catplot 'kind' input is allowed.
                Defaults to "box".
            what (str, optional): What to plot. For trajectories valid inputs are
                sigma or temperature. Otherwise, it's the keys of the stats dict.
                I.e., t_nucleation, T_nucleation, t_solidification.
                Defaults to "temperature".
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            context (bool, optional): Whether or not to print additional context
                in the figure. Mainly used to add shelf temperature.
        """
        stats_df, traj_df = self.to_frame()
        if kind.lower().startswith("traj"):
            if self._emptyStore:
                raise ValueError("No trajectories have been stored. Cannot plot.")
            df = traj_df
        else:
            df = stats_df

        if group != "all":
            if not isinstance(group, (tuple, list)):
                group = [group]
            df = df[df.group.isin(group)]

        if kind.lower().startswith("traj"):
            df = df[df.state.str.contains(what)]
            if len(df) == 0:
                raise ValueError(
                    f"No data for {what} for plot type {kind} and group {group}. "
                    + "Please check your inputs."
                )
            ax = sns.lineplot(data=df, hue="group", y="value", x="Time")

            if (what == "temperature") and context:
                N_timeSteps = (
                    int(np.ceil(self.opcond.t_tot / self.dt)) + 1
                )  # the total number of timesteps
                t = np.arange(0, N_timeSteps * self.dt, self.dt)
                T_shelf = self.opcond.tempProfile(self.dt)
                sns.lineplot(
                    x=t,
                    y=T_shelf,
                    ax=ax,
                    linestyle="dashed",
                    linewidth=0.75,
                    color="black",
                    label="Shelf Temperature",
                )

        else:
            df = df[df.variable.str.contains(what)]
            if len(df) == 0:
                raise ValueError(
                    f"No data for {what} for plot type {kind} and group {group}. "
                    + "Please check your inputs."
                )
            sns.catplot(data=df, hue="group", y="value", kind=kind, x="variable")

    def _sigmaCrossingIndices(
        self, threshold: float = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return indices where sigma crosses threshold.

        A helper function identifying the time (indices) in the state matrix
        where sigma crosses a certain value.
        Args:
            threshold (Optional[float], optional): The threshold to apply.
                Defaults to self.solidificationThreshold.

        Raises:
            ValueError: Simulation needs to be run first.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Indices where threshold is reached and
                mask identifying vials that never reach threshold.
        """
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run first.")

        if threshold is None:
            threshold = self.solidificationThreshold

        Indices = np.argmax(self.X_sigma > threshold, axis=1)
        # vials that never exceeded solidification threshold
        neverReached = ~np.any(self.X_sigma > threshold, axis=1)

        return Indices, neverReached

    def _buildInteractionMatrices(self) -> Tuple[csr_matrix, np.ndarray]:
        """Build interaction matrices.

        Returns:
            Tuple[csr_matrix, np.ndarray]: Internal interaction matrix and external
                interaction vector.
        """
        # This is a (n_x*n_y*n_z) x (n_x*n_y*n_z) matrix
        # of interactions between vial pairs
        # Please note that, on any given shelf, we index vials this way:
        # [[1, 2, 3, 4, ..., n_x],
        # [n_x+1, n_x+2, ..., 2*n_x],
        # [...],
        # [(n_y-1)*n_x+1, ... , n_y*n_x]]
        # If one thinks of a shelf as an array, we hence use 'row-major' ordering.
        # This is the same as numpy (most of the time), but different from Matlab.
        # These matrices can then be used as a sort of mask to overlay on, e.g.,
        # temperatures to calculate the correct driving forces for each vial/pair.

        n_x = self.N_vials[0]  # the number of vials in the horizontal direction
        n_y = self.N_vials[1]  # the number of vials in the vertical direction
        n_z = self.N_vials[2]  # the number of vials in the upwards/downwards direction
        # (perpendicular to the plane spanned by x and y)

        # create interaction matrix for horizontal (x-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,2), (2,3), ..., (i,i+1) have horizontal interactions,
        # except where i = n_x!
        # this means that in the IA matrix off-diagonal elements are 1
        # except where i%n_x==0 or j%n_x==0
        dx_pattern = np.ones((n_x * n_y * n_z - 1,))
        # the vial at index position i+1 is at the edge -> one neighbor less
        idx_delete = [i for i in range(n_x * n_y * n_z - 1) if (i + 1) % n_x == 0]
        dx_pattern[idx_delete] = 0
        # we store this as a compressed sparse row (most efficient format?)
        DX = csr_matrix(np.diag(dx_pattern, k=1) + np.diag(dx_pattern, k=-1))

        # create interaction matrix for vertical (y-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,n_x+1), (2,n_x+2), ..., (i,n_x+i) have vertical interactions
        # for i in [1, n_x*(n_y-1)]
        # this means that in the IA matrix n_x-removed-off-diagonal elements are 1
        dy_pattern = np.ones((n_x * (n_y * n_z - 1),))
        idy_delete = [
            i for i in range(n_x * (n_y * n_z - 1) - 1) if (i + 1) % (n_x * n_y) == 0
        ]
        for i in range(n_x):  # need to extract n_x points (i.e. an entire row)
            # every time we have an interaction
            delete_n_x = [k - i for k in idy_delete]  # remove all interactions
            dy_pattern[delete_n_x] = 0

        if self.const["vial_arrangement"] == "square":
            DY = csr_matrix(np.diag(dy_pattern, k=n_x) + np.diag(dy_pattern, k=-n_x))

        else:
            # for hexagonal arrangement matrix contains more elements, hence additional dy_pattern
            # elements, also we build the DY matrix for each z-layer and then construct it
            # later for vertical directions (AK, 24.04.2024)
            # COnfiguration details: (1) in z-direction, the vials are directly on top of each other, no change of positions
            # across layers, (2) the first vial in the top left corner (1,1) interact in y-direction with (2,1),
            # meaning that the second row is shifted left to create the hexagonal arrangement
            # TODO: refine the definition of dy_patterns and documentation

            dy_pattern_center = np.ones((n_x * (n_y - 1),))
            if n_x * (n_y - 1) > 0:
                dy_pattern_extra_upper = np.ones((n_x * (n_y - 1) - 1,))
                dy_pattern_extra_lower = np.ones((n_x * (n_y - 1) + 1,))

                if n_x * (n_y - 1) > 1:
                    dy_pattern_extra_upper[0] = 0
                dy_pattern_extra_lower[0] = 0

                for j in range(n_y):
                    if j % 2 == 0:
                        if j == 0:
                            dy_pattern_extra_upper[j * n_x + 1 : (j + 1) * n_x] = 0
                        else:
                            dy_pattern_extra_upper[j * n_x - 1 : (j + 1) * n_x] = 0

                for j in range(n_y):
                    if j % 2 == 1:
                        dy_pattern_extra_lower[j * n_x : (j + 1) * n_x + 1] = 0

                        # print(j%2, j * n_x,  (j + 1) * n_x + 1)

                # for i in range(n_x):  # need to extract n_x points (i.e. an entire row)
                #     # every time we have an interaction
                #     delete_n_x = [k - i for k in idy_delete]  # remove all interactions
                #     dy_pattern_extra_upper[delete_n_x] = 0
                #     dy_pattern_extra_lower[delete_n_x] = 0

                # print("\n dy_pattern_center: ", dy_pattern_center)
                # print("\n dy_pattern_extra_upper: ", dy_pattern_extra_upper)
                # print("\n dy_pattern_extra_lower: ", dy_pattern_extra_lower)

                DY_one_layer = (
                    np.diag(dy_pattern_center, k=n_x)
                    + np.diag(dy_pattern_center, k=-n_x)
                    + np.diag(dy_pattern_extra_upper, k=n_x + 1)
                    + np.diag(dy_pattern_extra_lower, k=n_x - 1)
                    + np.diag(dy_pattern_extra_upper, k=-(n_x + 1))
                    + np.diag(dy_pattern_extra_lower, k=-(n_x - 1))
                )

                # print("\n DY_one_layer: \n", DY_one_layer)

                DY_total = np.zeros((n_x * n_y * n_z, n_x * n_y * n_z))
                for i in range(n_z):
                    DY_total[
                        i * (n_x * n_y) : (i + 1) * (n_x * n_y),
                        i * (n_x * n_y) : (i + 1) * (n_x * n_y),
                    ] = DY_one_layer

                DY = csr_matrix(DY_total)

        # create interaction matrix for upwards/downwards direction interactions
        # matrix is 1 where two vials have an interaction and 0 otherwise
        # pairs (1,(n_x*n_y)+1), (2,(n_x+n_y)+2), ..., (i,i+(n_x+n_y)) have interactions
        # for i in [1, n_x*n_y*(n_z-1)]
        if n_z > 1:
            dz_pattern = np.ones(
                (n_x * n_y * (n_z - 1),)
            )  # @DRO: For n_z = 1, there will be no pattern
            DZ = csr_matrix(
                np.diag(dz_pattern, k=(n_x * n_y)) + np.diag(dz_pattern, k=(-n_x * n_y))
            )
        else:
            # no pattern in z-direction
            DZ = 0

        # how many interactions does each vial
        # have with other vials
        # diagflat because sum over a sparse matrix
        # returns a np.matrix object, not an ndarray
        VIAL_INT = csr_matrix(np.diagflat(np.sum(DX + DY + DZ, axis=1)))

        # the overall interaction matrix is given by the sum
        # of the interaction matrices DX and DY
        # minus the diagonal containing the sum of all vial-vial interactions
        interactionMatrix = DX + DY + DZ - VIAL_INT

        # at most, any cubic vial on a 2D shelf can have
        # 4 interactions. 4 - VIAL_INT is the
        # number of external interactions (excl. the shelf).
        # For the 3D case this max is 6.

        # max interactions depends on the vial arrangement
        if self.const["vial_arrangement"] == "square":
            maxInteractions = 4 + (n_z > 1) * 2
        else:
            maxInteractions = 6 + (n_z > 1) * 2

        VIAL_EXT = maxInteractions * np.ones((n_x * n_y * n_z,)) - VIAL_INT.diagonal()

        return interactionMatrix, VIAL_EXT

    def _buildShelfHeatFlow(
        self,
    ):
        """Build the shelf heat flow array.

        Because the shelf heat flow may be dependent on the rng if
        s_sigma_rel > 0 we need to separate this process out so we
        can reinitialize this particular array for Snowfall.
        """
        # store which seed was used to compute this
        # (getter method will know to recall this function
        # if the value changes)
        self._seedUsed = self.seed

        if self.N_vials[2] == 1:
            # 2D/shelf simulation
            if ("s_sigma_rel" in self.k.keys()) and (self.k["s_sigma_rel"] > 0):
                self.k["shelf"] = (
                    self.k["s0"]
                    + self._rng.normal(size=self.N_vials_total)
                    * self.k["s_sigma_rel"]
                    * self.k["s0"]
                )
            else:
                self.k["shelf"] = self.k["s0"]
        else:
            # 3D/pallet simulation
            self.k["shelf"] = 0

        if np.any(self.k["shelf"] < 0):
            # k cannot be smaller than 0
            self.k["shelf"][self.k["shelf"] < 0] = 0
            print("There were shelf heat transfer coefficients < 0. I set them to 0.")

        self._H_shelf = self.k["shelf"] * self.const["A"]  # either a scalar or a vector

    def _buildHeatflowMatrices(self):
        """Build heatflow matrices."""
        # store which N_vials was used to compute this
        # (getter method will know to recall this function
        # if the value changes)
        self._NvialsUsed = self.N_vials

        A = self.const["A"]

        interactionMatrix, VIAL_EXT = self._buildInteractionMatrices()

        # compute H_int
        self._H_int = interactionMatrix * self.k["int"] * A

        # compute H_ext
        self._H_ext = VIAL_EXT * self.k["ext"] * A

        # compute H_shelf (potentially depends on rng/seed)
        # will be zero in case of 3D system (assumes pallet in storage freezing)
        self._buildShelfHeatFlow()

    def to_frame(
        self, n_timeSteps: int = 250
    ) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        """Save statistics (and states) as pandas dataframe.

        Args:
            n_timeSteps (int, optional): For states a reduced representation is stored.
                This defines the number of subsamples. Defaults to 250.

        Raises:
            ValueError: Simulation needs to run first.

        Returns:
            Tuple[pd.DataFrame, Optional[pd.DataFrame]]: The statistics and if available
                the states over time (both in long form).
        """
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run first.")

        df = pd.DataFrame(self.stats)
        df.index.name = "vial"
        df = df.reset_index()

        _, VIAL_EXT = self._buildInteractionMatrices()

        df["group"] = VIAL_EXT
        # square arrangement
        if self.const["vial_arrangement"] == "square":
            df.loc[df.group == 3 - (self.N_vials[2] == 1), "group"] = "corner"
            df.loc[df.group == 2 - (self.N_vials[2] == 1), "group"] = "edge"
            df.loc[df.group == 1, "group"] = "side"
            df.loc[df.group == 0, "group"] = "core"
        # hexagonal arrangement
        else:
            df.loc[df.group == 5 - (self.N_vials[2] == 1), "group"] = "corner"
            df.loc[
                (df.group == 1)
                | (df.group == 2)
                | (df.group == 3)
                | (df.group == 4 - (self.N_vials[2] == 1)),
                "group",
            ] = "edge"
            df.loc[
                (df.group == 1)
                | (df.group == 2)
                | (df.group == 3)
                | (df.group == 4 - (self.N_vials[2] == 1)),
                "group",
            ] = "side"
            df.loc[df.group == 0, "group"] = "core"

        stats_df = df.melt(id_vars=["group", "vial"])

        n_storedStates = np.sum(self._storageMask)
        if n_storedStates > 0:
            # to reduce computational cost we limit number of timepoints to n_timeSteps
            df = pd.DataFrame(self._X[:, :: int(self._X.shape[1] / (n_timeSteps - 1))])
            df.columns = self._t[:: int(self._X.shape[1] / (n_timeSteps - 1))]
            df.columns.name = "Time"

            df["state"] = np.concatenate(
                [
                    np.repeat("temperature", n_storedStates),
                    np.repeat("sigma", n_storedStates),
                ]
            )
            df["vial"] = np.tile(np.where(self._storageMask)[0], 2)
            df["group"] = np.tile(VIAL_EXT[self._storageMask], 2)
            # square arrangement
            if self.const["vial_arrangement"] == "square":
                df.loc[df.group == 3 - (self.N_vials[2] == 1), "group"] = "corner"
                df.loc[df.group == 2 - (self.N_vials[2] == 1), "group"] = "edge"
                df.loc[df.group == 1 - (self.N_vials[2] == 1), "group"] = "side"
                df.loc[df.group == 0, "group"] = "core"

            # hexagonal arrangement
            else:
                df.loc[df.group == 5 - (self.N_vials[2] == 1), "group"] = "corner"
                df.loc[
                    (df.group == 1)
                    | (df.group == 2)
                    | (df.group == 3)
                    | (df.group == 4 - (self.N_vials[2] == 1)),
                    "group",
                ] = "edge"
                df.loc[
                    (df.group == 1)
                    | (df.group == 2)
                    | (df.group == 3)
                    | (df.group == 4 - (self.N_vials[2] == 1)),
                    "group",
                ] = "side"
                df.loc[df.group == 0, "group"] = "core"

            traj_df = df.melt(id_vars=["group", "vial", "state"])
        else:
            traj_df = None

        return stats_df, traj_df

    def __repr__(self) -> str:
        """Return string representation of the Snowflake class.

        A simple indicator of the most important properties of the Snowflake.
        Returns:
            str: The Snowflake class string representation giving some basic info.
        """
        return (
            f"Snowflake([N_vials: {self.N_vials}, "
            + f"dt: {self.dt}, seed: {self.seed}, seed_v: {self.seed_v}])"
        )
