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

        self.pool = mp.Pool(self.pool_size)
        self.Nrep = Nrep

        if 'seed' in kwargs.keys():
            del kwargs['seed']

        # fix seeds
        self.snowflakes = [Snowflake(seed=i, *kwargs) for i in range(int(Nrep))]

    def run(self):
        # run the individual snowflakes in a parallelized manner
        # processes = []
        # for sf in self.snowflakes:
        #     p = mp.Process(target=sf.run)
        #     processes.append(p)
        # print(processes)
        # [x.start() for x in processes]
        # print(self.snowflakes[0], 'main')
        [sf.run() for sf in self.snowflakes]

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
