'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

from ethz_snow.operatingConditions import OperatingConditions
from ethz_snow.constants import (
    A, hl, T_eq, kb, V, b, cp_solution,
    mass, solid_fraction, cp_s, depression,
    alpha, beta_solution, cp_i, cp_w)

import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.sparse import csr_matrix, lil, lil_matrix
from typing import List, Tuple, Union, Sequence, Optional

from ethz_snow.snowflake import Snowflake

HEATFLOW_REQUIREDKEYS = ('int', 'ext', 's0')
VIAL_GROUPS = ('corner', 'edge', 'center', 'all')


class Snowfall():
    def __init__(
        self,
        Nrep: int = 5e3,
        pool_size: int = 12,
        **kwargs
    ):

        self.pool_size = pool_size
        self.Nrep = Nrep

        # fix seeds
        self.snowflakes = [Snowflake(*kwargs) for _ in range(Nrep)]

    def run(self):
        # run the individual snowflakes in a parallelized manner

    def nucleationTimes(self, group: Union[str, Sequence[str]] = 'all',
                        fromStates: bool = False) -> np.ndarray:

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
            t_nucleation = self.stats['t_nucleation']

        I_groups = self.getVialGroup(group)

        return t_nucleation[I_groups]

    def nucleationTemperatures(self, group: Union[str, Sequence[str]] = 'all',
                               fromStates: bool = False) -> np.ndarray:

        if fromStates:
            I_nucleation, lateBloomers = self._sigmaCrossingIndices(threshold=0)
            T_nucleation_states = (np.array([self.X_T[i, I-1]
                                   for i, I in zip(range(self.N_vials_total), I_nucleation)])
                                   .astype(float))  # should always be float, but just to be sure

            # non-nucleated vials are set to NaN
            T_nucleation_states[lateBloomers] = np.nan

            T_nucleation = np.full(self.N_vials_total, np.nan)
            T_nucleation[self._storageMask] = T_nucleation_states

        else:
            T_nucleation = self.stats['T_nucleation']

        I_groups = self.getVialGroup(group)

        return T_nucleation[I_groups]

    def solidificationTimes(self, group: Union[str, Sequence[str]] = 'all', threshold: float = 0.9,
                            fromStates: bool = False) -> np.ndarray:

        if fromStates:
            t_nucleation = self.nucleationTimes(fromStates=fromStates)

            I_solidification, neverGrownUp = self._sigmaCrossingIndices(threshold=threshold)

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
            t_solidification = self.stats['t_solidification']

        I_groups = self.getVialGroup(group)

        return t_solidification[I_groups]

    def sigmaCounter(self,
                     time: Union[Sequence[float], float],
                     threshold: float = 0.9,
                     fromStates: bool = False) -> np.ndarray:
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run before induction times can be extracted.")

        if isinstance(time, (float, int)):
            time = [time]

        counter = np.zeros(len(time))
        for i, t in enumerate(time):

            if fromStates:
                I_time = np.argmax(self._t >= t)
                counter[i] = np.sum(self.X_sigma[:, I_time] > threshold, axis=0)
            else:
                counter[i] = np.sum(self.stats['t_solidification'] <= t)

        return counter

    def plot(self,
             what: str = 'temperature', kind: str = 'trajectories',
             group: Union[str, Sequence[str]] = 'all'):

        stats_df, traj_df = self.toDataframe()
        if not kind.lower().startswith('traj'):
            df = stats_df
        else:
            df = traj_df

        if group != 'all':
            if not isinstance(group, (tuple, list)):
                group = [group]
            df = df[df.group.isin(group)]

        if not kind.lower().startswith('traj'):
            df = df[df.variable.str.lower().str.contains(what.lower())]
            sns.catplot(data=df, hue='group', y='value', kind=kind, x='variable')
        else:
            df = df[df.state.str.lower().str.contains(what.lower())]
            sns.lineplot(data=df, hue='group', y='value', x='Time')

    def toDataframe(self, n_timeSteps=250) -> Tuple[np.ndarray, Optional[np.ndarray]]:

        # does it make sense to have this?

    def __repr__(self) -> str:
        """ The string representation of the Snowfall class.

        Returns:
            str: The Snowfall class string representation giving some basic info.
        """

        return (f"Snowfall([{self.Nrep} Snowflakes, pool_size: {self.pool_size}])")
