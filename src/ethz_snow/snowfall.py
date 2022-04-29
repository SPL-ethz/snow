"""Implement Snowfall class.

This module contains the Snowfall class used to run repeated (!)
simulations of water nucleation in vials. It makes use of class
Snowflake for the individual simulations.
"""
import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns

import multiprocessing as mp
from typing import List, Tuple, Union, Sequence, Optional

from ethz_snow.snowflake import Snowflake


class Snowfall:
    """A class to handle multiple Stochastic Nucleation of Water simulation.

    More information regarding the equations and their derivation can be found in
    "Stochastic shelf-scale modeling framework for the freezing stage in
    freeze-drying processes", Deck, Ochsenbein, and Mazzotti (2022),
    Int J Pharm, 613, 121276, https://doi.org/10.1016/j.ijpharm.2021.121276.

    Parameters:
        Nrep (int): Number of repetitions.
        pool_size (int): Size of worker pool for parallelization.
        stats (dict): Statistics for each simulation.
        stats_df (pd.DataFrame): Long-form table of all statistics.
        simulationStatus (int): Status of simulation (0 = not run, 1 = run).
    """

    def __init__(self, Nrep: int = 5, pool_size: int = None, **kwargs):
        """Construct a Snowfall object.

        Args:
            Nrep (int, optional): Number of repetitions. Defaults to 5.
            pool_size (int, optional): Size of worker pool for parallelization.
                Defaults to None (= # of cpu cores).
        """
        self.pool_size = pool_size

        self.Nrep = int(Nrep)

        if "seed" in kwargs.keys():
            # seed will be chosen by Snowfall
            del kwargs["seed"]

        if ("storeStates" in kwargs.keys()) and (kwargs["storeStates"] is not None):
            print("WARNING: States cannot be stored for Snowfall simulations.")
            del kwargs["storeStates"]

        Sf_template = Snowflake(**kwargs)
        Sf_template._buildHeatflowMatrices()  # pre-build H_int, H_ext, H_shelf

        self.Sf_template = Sf_template

        self.stats = dict()
        self.stats_df = pd.DataFrame()

    @property
    def simulationStatus(self) -> int:
        """Return simulation status of instance.

        Returns:
            int: Simulation status. 0 = not run, 1 = run.
        """
        if self.stats:
            return 1
        else:
            return 0

    @classmethod
    def _uniqueFlake(cls, S, seed):
        """Run a single Snowflake with a specific seed."""
        S.seed = seed
        S.run()
        return S.stats

    @classmethod
    def _uniqueFlake_sync(cls, S, seed, return_dict):
        """Run a single Snowflake with a specific seed.

        Run a single Snowflake with a specific seed
        and write it back into multiprocessing.Manager object for data sharing.
        """
        S.seed = seed
        S.run()
        return_dict[seed] = S.stats

    def run(self, how="async"):
        """Run all the Snowflake simulations.

        Args:
            how (str, optional): How to perform runs. Valid options are
                'async', 'sync', and 'sequential' (no parallelization).
                Defaults to "async".
        """
        # clean up old simulation
        self.stats = dict()
        self.stats_df = pd.DataFrame()

        if how == "async":
            with mp.Pool(self.pool_size) as p:
                # starmap is only available since python 3.3
                # it allows passing multiple arguments
                res = p.starmap_async(
                    Snowfall._uniqueFlake,
                    [(self.Sf_template, i) for i in range(self.Nrep)],
                ).get()
            self.stats = {i: r for i, r in enumerate(res)}

        elif how == "sync":
            manager = mp.Manager()
            return_dict = manager.dict()
            with mp.Pool(self.pool_size) as p:
                # starmap is only available since python 3.3
                # it allows passing multiple arguments
                p.starmap(
                    Snowfall._uniqueFlake_sync,
                    [(self.Sf_template, i, return_dict) for i in range(self.Nrep)],
                )
            self.stats = dict(return_dict)

        elif how == "sequential":
            for i in range(self.Nrep):
                self.stats[i] = self._uniqueFlake(self.Sf_template, i)

    def nucleationTimes(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:
        """Return nucleation times.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            seed (Union[int, Sequence[int], None], optional): Seed(s) to return.
                Defaults to None.

        Returns:
            np.ndarray: The nucleation times.
        """
        return self._returnStats(what="tnuc", group=group, seed=seed)

    def nucleationTemperatures(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:
        """Return nucleation temperatures.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            seed (Union[int, Sequence[int], None], optional): Seed(s) to return.
                Defaults to None.

        Returns:
            np.ndarray: The nucleation temperatures.
        """
        return self._returnStats(what="Tnuc", group=group, seed=seed)

    def solidificationTimes(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:
        """Return solidification times.

        Args:
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            seed (Union[int, Sequence[int], None], optional): Seed(s) to return.
                Defaults to None.

        Returns:
            np.ndarray: The solidification times.
        """
        return self._returnStats(what="tsol", group=group, seed=seed)

    def _returnStats(
        self,
        what: str,
        group: str = "all",
        seed: Union[int, Sequence[int], None] = None,
    ):
        """Return arbitarary statistic.

        Helper function to compute arbitrary statistic for
        nucleationTimes, solidificationTimes, nucleationTemperatures.
        Args:
            what (str): What to return.
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
            seed (Union[int, Sequence[int], None], optional): Seed(s) to return.
                Defaults to None.

            Returns:
                np.ndarray: The statistic in question.
        """
        df = self.to_frame()

        if group != "all":
            if not isinstance(group, (tuple, list)):
                group = [group]
            df = df[df.group.isin(group)]

        if seed is not None:
            if not isinstance(seed, (tuple, list)):
                seed = [seed]
            df = df[df.seed.isin(seed)]

        if what == "tnuc":
            df = df[df.variable == "t_nucleation"]
        elif what == "Tnuc":
            df = df[df.variable == "T_nucleation"]
        elif what == "tsol":
            df = df[df.variable == "t_solidification"]

        return df["value"].to_numpy()

    def plot(
        self,
        what: str = "t_nucleation",
        kind: str = "box",
        seed: Union[int, Sequence[int], None] = None,
        group: Union[str, Sequence[str]] = "all",
    ):
        """Create plots for Snowfall object.

        Args:
            kind (str, optional): Any sns.catplot 'kind' input is allowed.
                Defaults to "box".
            what (str, optional): What to plot, i.e., keys of the stats dict.
                Valid options are t_nucleation, T_nucleation, t_solidification.
                Defaults to "t_nucleation".
            seed (Union[int, Sequence[int], None], optional): Seed(s) to return.
                Defaults to None.
            group (Union[str, Sequence[str]], optional): Subgroup to return.
                Defaults to "all".
        """
        df = self.to_frame()

        if group != "all":
            if not isinstance(group, (tuple, list)):
                group = [group]
            df = df[df.group.isin(group)]

        if seed is not None:
            if not isinstance(seed, (tuple, list)):
                seed = [seed]
            df = df[df.seed.isin(seed)]

        df = df[df.variable.str.contains(what)]
        sns.catplot(data=df, hue="group", y="value", kind=kind, x="variable")

    def to_frame(self) -> pd.DataFrame:
        """Save statistics as pandas dataframe.

        Raises:
            ValueError: Simulation needs to run first.

        Returns:
            pd.DataFrame: The statistics in long form.
        """
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run first.")
        elif self.stats_df.empty:
            stats_df = pd.DataFrame(
                columns=["group", "vial", "variable", "value", "seed"]
            )
            self.Sf_template._simulationStatus = 1
            for i in range(self.Nrep):
                self.Sf_template.stats = self.stats[i]
                loc_stats_df, _ = self.Sf_template.to_frame()
                loc_stats_df["seed"] = i
                stats_df = stats_df.append(loc_stats_df)

            # clean up template
            self.Sf_template.stats = dict()
            self.Sf_template._simulationStatus = 0

            self.stats_df = stats_df
        elif not self.stats_df.empty:
            # do not recalculate, but load from memory
            stats_df = self.stats_df

        return stats_df

    def __repr__(self) -> str:
        """Return string representation of the Snowfall class.

        Returns:
            str: The Snowfall class string representation giving some basic info.
        """
        return (
            f"Snowfall([{self.Nrep} Snowflake{'s' if self.Nrep > 1 else ''}, "
            + f"pool_size: {self.pool_size}])"
        )
