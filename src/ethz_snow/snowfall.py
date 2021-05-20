import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns

import multiprocessing as mp
from typing import List, Tuple, Union, Sequence, Optional

from ethz_snow.snowflake import Snowflake


class Snowfall:
    def __init__(self, Nrep: int = 5, pool_size: int = None, **kwargs):

        self.pool_size = pool_size

        # self.pool = mp.Pool(self.pool_size)
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
    def simulationStatus(self):
        if self.stats:
            return 1
        else:
            return 0

    @classmethod
    def uniqueFlake(cls, S, seed):
        S.seed = seed
        S.run()
        return S.stats

    @classmethod
    def uniqueFlake_sync(cls, S, seed, return_dict):
        S.seed = seed
        S.run()
        return_dict[seed] = S.stats

    def run(self, how="async"):

        # clean up old simulation
        self.stats = dict()
        self.stats_df = pd.DataFrame()

        # run the individual snowflakes in a parallelized manner
        if how == "async":
            with mp.Pool(self.pool_size) as p:
                # starmap is only available since python 3.3
                # it allows passing multiple arguments
                res = p.starmap_async(
                    Snowfall.uniqueFlake,
                    [(self.Sf_template, i) for i in range(self.Nrep)],
                ).get()
            self.stats = res

        elif how == "sync":
            manager = mp.Manager()
            return_dict = manager.dict()
            with mp.Pool(self.pool_size) as p:
                # starmap is only available since python 3.3
                # it allows passing multiple arguments
                p.starmap(
                    Snowfall.uniqueFlake_sync,
                    [(self.Sf_template, i, return_dict) for i in range(self.Nrep)],
                )
            self.stats = dict(return_dict)

        elif how == "sequential":
            for i in range(self.Nrep):
                self.stats[i] = self.uniqueFlake(self.Sf_template, i)

    def nucleationTimes(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:
        return self._returnStats(what="tnuc", group=group, seed=seed)

    def nucleationTemperatures(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:
        return self._returnStats(what="Tnuc", group=group, seed=seed)

    def solidificationTimes(
        self,
        group: Union[str, Sequence[str]] = "all",
        seed: Union[int, Sequence[int], None] = None,
    ) -> np.ndarray:

        return self._returnStats(what="tsol", group=group, seed=seed)

    def _returnStats(
        self,
        what: str,
        group: str = "all",
        seed: Union[int, Sequence[int], None] = None,
    ):
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
            df = df[df.variable == 't"solidification']

        return df["value"].to_numpy()

    def plot(
        self,
        what: str = "t_nucleation",
        kind: str = "box",
        seed: Union[int, Sequence[int], None] = None,
        group: Union[str, Sequence[str]] = "all",
    ):
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
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run first.")
        elif self.stats_df.empty:
            stats_df = pd.DataFrame(
                columns=["group", "vial", "variable", "value", "seed"]
            )
            for i in range(self.Nrep):
                self.Sf_template.stats = self.stats[i]
                loc_stats_df, _ = self.Sf_template.to_frame()
                loc_stats_df["seed"] = i
                stats_df = stats_df.append(loc_stats_df)

            # clean up template
            self.Sf_template.stats = dict()

            self.stats_df = stats_df

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
