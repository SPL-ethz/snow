"""Implement Snowflake class.

This module contains the Snowflake class used to run simulations
of water nucleation in vials.
"""
from ethz_snow.operatingConditions import OperatingConditions
from ethz_snow.constants import (
    A,
    hl,
    T_eq,
    kb,
    V,
    b,
    cp_solution,
    mass,
    solid_fraction,
    cp_s,
    depression,
    alpha,
    beta_solution,
    cp_i,
    cp_w,
)

import numpy as np
import pandas as pd
import re
import sys

# import matplotlib.pyplot as plt
import seaborn as sns

from scipy.sparse import csr_matrix
from typing import List, Tuple, Union, Sequence, Optional

import time

HEATFLOW_REQUIREDKEYS = ("int", "ext", "s0")
VIAL_GROUPS = ("corner", "edge", "center", "all")


class Snowflake:
    """A class to handle a single Stochastic Nucleation of Water simulation.

    More information regarding the equations and their derivation can be found in
    XXX, Deck et al. (2021).

    Attributes:
        H_ext (np.ndarray): External heat transfer vector.
        H_int (csr_matrix): Internal heat transfer matrix.
        H_shelf (np.ndarray): Shelf heat transfer vector.
        N_vials (tuple): Number of vials in each dimension.
        N_vials_total (int): Total number of vials.
        X_T (np.ndarray): Temperature state over time.
        X_sigma (np.ndarray): Sigma state over time.
        dt (float): Time step.
        k (dict): Heat transfer coefficients.
        seed (int): Seed to be used in rng.
        solidificationThreshold (float): What sigma value constitutes a 'solid'.
        stats (dict): Run statistics (nucleation time, etc.).

    """

    def __init__(
        self,
        k: dict = {"int": 20, "ext": 20, "s0": 20, "s_sigma_rel": 0.1},
        N_vials: Tuple[int, int, int] = (
            7,
            7,
            1,
        ),  # should this be part of operating conditions?
        storeStates: Union[str, Sequence[str], Sequence[int], None] = None,
        solidificationThreshold: float = 0.9,
        dt: float = 2,
        seed: int = 2021,
        opcond: OperatingConditions = OperatingConditions(),
    ):
        """Construct a Snowflake object.

        Args:
            k (dict, optional): A dictionary containing the heat transfer coefficients.
                Must contain keys 'int', 'ext', 's0'.
                Defaults to {"int": 20, "ext": 20, "s0": 20, "s_sigma_rel": 0.1}.
            N_vials (Tuple[int, int, int], optional): Number of vials in each dimension.
                Defaults to ( 7, 7, 1).
            storeStates (Union[str, Sequence[str], Sequence[int], None]):
                Whether and for which vials to store states (temp, sigma).
            solidificationThreshold (float, optional): What sigma value is
                a 'solid'. Defaults to 0.9.
            dt (float, optional): Time step size. Defaults to 2.
            seed (int, optional): Seed to use for rng. Defaults to 2021.
            opcond (OperatingConditions, optional): Operating conditions to apply.
                Defaults to OperatingConditions().

        Raises:
            TypeError: If k is not a dict.
            ValueError: If k does not contain all necessary keys.
            NotImplementedError: If N_vials[2] is not 1. Only 2D model is implemented.
            TypeError: If opcond is not of type operatingConditions.
            ValueError: If storeStates is not meaningful.

        Examples:
            >>> Sf = Snowflake()
            >>> Sf = Snowflake(Nvials=(4, 3, 1), dt = 5)
            >>> Sf = Snowflake(storeStates = ['edge_random_4', 'uniform.center.5'])
            >>> Sf = Snowflake(storeStates = (0, 4, 10))
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

        if isinstance(N_vials, list):
            # we exepct this to be immutable
            N_vials = tuple(N_vials)

        if N_vials[2] > 1:
            raise NotImplementedError(
                "Only 2D (shelf) models are implemented at this moment."
            )
        else:
            self.N_vials = N_vials

        self.dt = dt

        # store the seed to look it up if need be
        self.seed = seed

        if not isinstance(opcond, OperatingConditions):
            raise TypeError(
                "Input opcond must be an instance of class OperatingConditions."
            )
        else:
            self.opcond = opcond

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
                storageMask[storeStates] = True
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

        self.solidificationThreshold = solidificationThreshold
        self._X = None
        self._t = None
        self.stats = dict()

        self._simulationStatus = 0

        self._H_int = None
        self._H_ext = None
        self._H_shelf = None
        # remember what N_vials was used to build heat flow matrices
        # so if it changes we know to rebuild them
        self._NvialsUsed = self.N_vials

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
        if (self._H_shelf is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
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
        H_int = self.H_int
        H_ext = self.H_ext
        H_shelf = self.H_shelf

        N_timeSteps = (
            int(np.ceil(self.opcond.t_tot / self.dt)) + 1
        )  # the total number of timesteps
        N_vials_total = self.N_vials_total

        # Initial state of the system
        # initial temperature
        T_k = np.ones(N_vials_total) * self.opcond.cooling["start"]  # [C]
        # state of the vials
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

        k_CN = np.argmax(t >= self.opcond.cnt)

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
                # DRO: Leif, I've moved the dt out of the Hs
                # because I associate Q with fluxes
                T_k[liquidMask] = T_k[liquidMask] + q_k[liquidMask] / hl * self.dt
                # a mask that is True where the temperature is
                # below the equilibrium temperature
                superCooledMask = T_k < T_eq
                # a vial that is both liquid and supercooled can nucleate
                nucleationCandidatesMask = liquidMask & superCooledMask
                # the total number of nucleation candidates
                n_nucleationCandidates = np.sum(nucleationCandidatesMask)
                # toc3_l = time.perf_counter()
                # print(f"Time steps: {toc3_l-toc2b_l:4.2e}")

                # nucleation probabilities for candidates
                P = np.zeros(N_vials_total)
                diceRolls = np.zeros(N_vials_total)

                P[nucleationCandidatesMask] = (
                    kb * V * (T_eq - T_k[nucleationCandidatesMask]) ** b * self.dt
                )
                # when we reach the timepoint of controlled nucleation
                # all vials (that thermodynamically can) nucleate
                if k == k_CN:
                    P.fill(1)
                diceRolls[nucleationCandidatesMask] = self._rng.random(
                    n_nucleationCandidates
                )

                # Nucleation
                nucleatedVialsMask = nucleationCandidatesMask & (diceRolls < P)

                stats["t_nucleation"][nucleatedVialsMask] = t[k] + self.dt
                stats["T_nucleation"][nucleatedVialsMask] = T_k[nucleatedVialsMask]

                q0 = (T_eq - T_k[nucleatedVialsMask]) * cp_solution * mass

                sigma_k[nucleatedVialsMask] = -q0 / (alpha - beta_solution)
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
            raise ValueError(
                "Simulation needs to be run before induction times can be extracted."
            )

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

    def getVialGroup(self, group: Union[str, Sequence[str]] = "all") -> np.ndarray:
        """Return mask for given group.

        mask[i] = True iff vial[i] is in G where G can be a list of groups.
        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".

        Raises:
            ValueError: Group is not known.

        Returns:
            np.ndarray: A boolean mask of size (N_vials_tot,).

        Examples:
            >>> S = Snowflake(N_vials = (2, 2, 1))
            >>> S.getVialGroup('corner')
            array([True, True, True, True])

            >>> S.getVialGroup(['edge', 'center'])
            array([False, False, False, False])
        """
        if isinstance(group, str):
            group = [group]

        _, VIAL_EXT = self._buildInteractionMatrices()

        myMask = np.zeros(self.N_vials_total, dtype=bool)
        for g in group:
            if g == "corner":
                myMask = myMask | (VIAL_EXT == 2)
            elif g == "edge":
                myMask = myMask | (VIAL_EXT == 1)
            elif g == "center":
                myMask = myMask | (VIAL_EXT == 0)
            elif g == "all":
                myMask = np.ones(self.N_vials_total, dtype=bool)
                break
            else:
                raise ValueError(f"Group {g} is not a known vial group.")

        return myMask

    def plot(
        self,
        kind: str = "trajectories",
        what: str = "temperature",
        group: Union[str, Sequence[str]] = "all",
    ):
        """Create plots for Snowflake object.

        Args:
            kind (str, optional): Type of plot to print. If input is
                'trajectories' will print temperature or sigma profile.
                Otherwise will print from stats dict. In that case,
                any sns.catplot 'kind' input is allowed.
                Defaults to "trajectories".
            what (str, optional): What to plot. For trajectories valid inputs are
                sigma or temperature. Otherwise, it's the keys of the stats dict.
                I.e., t_nucleation, T_nucleation, t_solidification.
                Defaults to "temperature".
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
        """
        stats_df, traj_df = self.to_frame()
        if not kind.lower().startswith("traj"):
            df = stats_df
        else:
            df = traj_df

        if group != "all":
            if not isinstance(group, (tuple, list)):
                group = [group]
            df = df[df.group.isin(group)]

        if not kind.lower().startswith("traj"):
            df = df[df.variable.str.contains(what)]
            sns.catplot(data=df, hue="group", y="value", kind=kind, x="variable")
        else:
            df = df[df.state.str.contains(what)]
            sns.lineplot(data=df, hue="group", y="value", x="Time")

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
            raise ValueError(
                "Simulation needs to be run before induction times can be extracted."
            )

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
        dx_pattern = np.ones((n_x * n_y - 1,))
        # the vial at index position i+1 is at the edge -> one neighbor less
        idx_delete = [i for i in range(n_x * n_y - 1) if (i + 1) % n_x == 0]
        dx_pattern[idx_delete] = 0
        # we store this as a compressed sparse row (most efficient format?)
        DX = csr_matrix(np.diag(dx_pattern, k=1) + np.diag(dx_pattern, k=-1))

        # create interaction matrix for vertical (y-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,n_x+1), (2,n_x+2), ..., (i,n_x+i) have vertical interactions
        # for i in [1, n_x*(n_y-1)]
        # this means that in the IA matrix n_x-removed-off-diagonal elements are 1
        dy_pattern = np.ones((n_x * n_y - n_x,))
        DY = csr_matrix(np.diag(dy_pattern, k=n_x) + np.diag(dy_pattern, k=-n_x))

        # how many interactions does each vial
        # have with other vials
        # diagflat because sum over a sparse matrix
        # returns a np.matrix object, not an ndarray
        VIAL_INT = csr_matrix(np.diagflat(np.sum(DX + DY, axis=1)))

        # the overall interaction matrix is given by the sum
        # of the interaction matrices DX and DY
        # minus the diagonal containing the sum of all vial-vial interactions
        interactionMatrix = DX + DY - VIAL_INT

        # at most, any cubic vial on a 2D shelf can have
        # 4 interactions. 4 - VIAL_INT is the
        # number of external interactions (excl. the shelf)
        # is it worth storing this as a sparse matrix? - DRO XXX
        VIAL_EXT = 4 * np.ones((n_x * n_y,)) - VIAL_INT.diagonal()

        return interactionMatrix, VIAL_EXT

    def _buildHeatflowMatrices(self):
        """Build heatflow matrices."""
        self._NvialsUsed = self.N_vials

        interactionMatrix, VIAL_EXT = self._buildInteractionMatrices()

        self._H_int = interactionMatrix * self.k["int"] * A

        self._H_ext = VIAL_EXT * self.k["ext"] * A

        if "s_sigma_rel" in self.k.keys():
            self.k["shelf"] = self.k["s0"] + self._rng.normal(size=self.N_vials_total)
        else:
            self.k["shelf"] = self.k["s0"]

        self._H_shelf = self.k["shelf"] * A  # either a scalar or a vector

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
            raise ValueError(
                "Simulation needs to be run before induction times can be extracted."
            )

        df = pd.DataFrame(self.stats)
        df.index.name = "vial"
        df = df.reset_index()

        _, VIAL_EXT = self._buildInteractionMatrices()

        df["group"] = VIAL_EXT
        df.loc[df.group == 2, "group"] = "corner"
        df.loc[df.group == 1, "group"] = "edge"
        df.loc[df.group == 0, "group"] = "center"

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

            df.loc[df.group == 2, "group"] = "corner"
            df.loc[df.group == 1, "group"] = "edge"
            df.loc[df.group == 0, "group"] = "center"

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
            + f"dt: {self.dt}, seed: {self.seed}])"
        )
