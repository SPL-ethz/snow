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

import multiprocessing as mp
from typing import List, Tuple, Union, Sequence, Optional

from ethz_snow.snowflake import Snowflake


class Snowfall():
    def __init__(
        self,
        Nrep: int = 5,
        pool_size: int = None,
        **kwargs
    ):

        self.pool_size = pool_size

        # self.pool = mp.Pool(self.pool_size)
        self.Nrep = int(Nrep)

        if 'seed' in kwargs.keys():
            # seed will be chosen by Snowfall
            del kwargs['seed']

        self.sf_kwargs = kwargs

    @classmethod
    def uniqueFlake(cls, sf_kwargs, seed, return_dict):
        S = Snowflake(*sf_kwargs, seed=seed)
        S.run()
        return_dict[seed] = S.stats

    def run(self):
        # run the individual snowflakes in a parallelized manner
        manager = mp.Manager()
        return_dict = manager.dict()
        with mp.Pool() as p:
            # starmap is only available since python 3.3
            # it allows passing multiple arguments
            p.starmap(Snowfall.uniqueFlake,
                      [(self.sf_kwargs, i, return_dict) for i in range(self.Nrep)])

        self.stats = dict(return_dict)

    def nucleationTimes(self, group: Union[str, Sequence[str]] = 'all',
                        fromStates: bool = False) -> np.ndarray:
        pass
        # XXX

    def nucleationTemperatures(self, group: Union[str, Sequence[str]] = 'all',
                               fromStates: bool = False) -> np.ndarray:
        pass
        # XXX

    def solidificationTimes(self, group: Union[str, Sequence[str]] = 'all', threshold: float = 0.9,
                            fromStates: bool = False) -> np.ndarray:
        pass
        # XXX

    def plot(self,
             what: str = 'temperature', kind: str = 'trajectories',
             group: Union[str, Sequence[str]] = 'all'):
        pass
        # XXX

    def toDataframe(self, n_timeSteps=250) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        pass
        # does it make sense to have this?

    def __repr__(self) -> str:
        """ The string representation of the Snowfall class.

        Returns:
            str: The Snowfall class string representation giving some basic info.
        """

        return (f"Snowfall([{self.Nrep} Snowflakes, pool_size: {self.pool_size}])")
