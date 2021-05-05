'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

from typing import List
from operatingConditions import OperatingConditions

class Snowfall:
    def __init__(
        k: dict = {'int':20, 'ext': 20, 's0':20, 's_sigma_rel':0.1},
        N_vials: List[int] = [7,7,1], # should this be part of operating conditions?
        Nrep: int = 5e3,
        dt: float = 2,
        pool_size: int = 12,
        opcond: OperatingConditions = OperatingConditions()
    ):
