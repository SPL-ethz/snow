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
import re, sys

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.sparse import csr_matrix, lil, lil_matrix
from typing import List, Tuple, Union, Sequence, Optional

HEATFLOW_REQUIREDKEYS = ('int', 'ext', 's0')
VIAL_GROUPS = ('corner', 'edge', 'center','all')


class Snowfall:
    def __init__(
        self,
        k: dict = {'int': 20, 'ext': 20, 's0': 20, 's_sigma_rel': 0.1},
        N_vials: Tuple[int, int, int] = (7, 7, 1),  # should this be part of operating conditions?
        storeStates: Optional[Union[str, Sequence[str], Sequence[int]]] = None,
        solidificationThreshold: float = 0.9,
        Nrep: int = 5e3,
        dt: float = 2,
        pool_size: int = 12,
        opcond: OperatingConditions = OperatingConditions()
    ):

        if not isinstance(k, dict):
            raise TypeError(f"Input k must be of type dict. Was given {type(k)}.")
        elif not all([key in k.keys() for key in HEATFLOW_REQUIREDKEYS]):
            raise ValueError((f"A required key was missing from dictionary k, specifically "
                             + f"{set(HEATFLOW_REQUIREDKEYS) - set(k.keys())}."))
        else:
            self.k = k

        if N_vials[2] > 1:
            raise NotImplementedError("Only 2D (shelf) models are implemented at this moment.")
        else:
            self.N_vials = N_vials

        self.Nrep = Nrep

        self.dt = dt

        self.pool_size = pool_size

        self.opcond = opcond

        if isinstance(storeStates, (list, tuple)):
            if all([isinstance(x, int) for x in storeStates]):
                storageMask = np.zeros(self.N_vials_total, dtype=bool)
                storageMask[storeStates] = True
            elif all([isinstance(x, str) for x in storeStates]):
                storageMasks = list(map(lambda x: self._interpretStorageString(x), storeStates))
                storageMask = np.logical_or.reduce(storageMasks)
        elif isinstance(storeStates, str):
            storageMask = self._interpretStorageString(storeStates)
        elif storeStates is None:
            storageMask = np.zeros(self.N_vials_total, dtype=bool)

        self._storageMask = storageMask

        self.solidificationThreshold = solidificationThreshold
        self._X = None
        self._t = None

        self._simulationStatus = 0

    @property
    def simulationStatus(self):
        if (self._simulationStatus == 1) or (self._X is not None):
            self._simulationStatus = 1

        return self._simulationStatus

    @property
    def N_vials_total(self):
        return int(np.prod(self.N_vials))

    @property
    def X_T(self):
        return self._X[:int(np.sum(self._storageMask)), :]

    @property
    def X_sigma(self):
        return self._X[int(np.sum(self._storageMask)):, :]

    def _interpretStorageString(self, myString):
        myString = myString.lower()
        if not any([word in myString for word in list(VIAL_GROUPS)+['random', 'uniform']]):
            raise ValueError("No valid vial group or key word used in storeStates.")
        elif any([word in myString for word in list(VIAL_GROUPS)]):
            for group in VIAL_GROUPS:
                if group in myString:
                    storageMask = self.getVialGroup(group)
                    break
        elif ('random' in myString) or ('uniform' in myString):
            # assume vial group 'all' is implied
            storageMask = np.ones(self.N_vials_total, dtype=bool)

        if ('random' in myString) or ('uniform' in myString):
            howMany = re.findall(r'\d+', myString)
            if len(howMany) == 0:
                # unclear how many to pick, default is 10%
                howMany = int(np.ceil(0.1*self.N_vials_total))
            elif len(howMany) == 1:
                howMany = int(howMany[0])
            elif len(howMany) > 1:
                raise ValueError("storeStates strings must contain one number at most.")

            I_candidates = np.where(storageMask)[0]
            if 'random' in myString:
                I_toStore = np.random.choice(I_candidates, size=howMany, replace=False)
            elif 'uniform' in myString:
                stepSize = int(np.ceil(len(I_candidates)/howMany))
                I_fromCandidates = np.arange(0, len(I_candidates), stepSize, dtype=int)
                I_toStore = I_candidates[I_fromCandidates]
            storageMask = np.zeros(self.N_vials_total, dtype=bool)
            storageMask[I_toStore] = True

        return storageMask

    def run(self):
        # clean up any potential old simulations
        self._X = None
        self._t = None

        N_timeSteps = int(np.ceil(self.opcond.t_tot / self.dt))+1  # the total number of timesteps

        # 1. build interaction matrices
        self._buildHeatflowMatrices()

        # 2. Obtain external temperature profile over time
        T_shelf = self.opcond.tempProfile(self.dt)
        T_ext = T_shelf  # need to make a switch so that this can be decoupled - DRO XX

        # 3. Initial state of the system
        # initial temperature
        T_k = np.ones((self.N_vials_total,)) * self.opcond.cooling['start']  # [C]
        # state of the vials
        sigma_k = np.zeros(self.N_vials_total)  # (0 = liq, 1 = completely frozen)

        # 4. Pre-allocate memory
        # our state matrix containing the entire history of the system
        X = np.zeros((2*np.sum(self._storageMask), N_timeSteps))
        t = np.arange(0, N_timeSteps * self.dt, self.dt)
        stats = dict(
            nucleationTimes=np.full(self.N_vials_total, np.nan),
            nucleationTemperatures=np.full(self.N_vials_total, np.nan),
            solidificationTimes=np.full(self.N_vials_total, np.nan))

        # 5. Iterate over time steps
        for k in np.arange(N_timeSteps):

            T_shelf_k = T_shelf[k]
            T_ext_k = T_ext[k]

            x_k = np.concatenate([T_k, sigma_k])
            X[:, k] = x_k[np.tile(self._storageMask, 2)]

            # a mask that is True where a vial is still fully liquid
            liquidMask = (sigma_k == 0)
            solidMask = ~liquidMask  # for convenience
            solidifiedMask = ((sigma_k > self.solidificationThreshold)
                              & np.isnan(stats['solidificationTimes']))
            stats['solidificationTimes'][solidifiedMask] = (t[k]
                                                            - stats['nucleationTimes'][solidifiedMask])

            # calculate heatflows for all vials
            q_k = (self.H_int @ T_k
                   + self.H_ext * (T_ext_k - T_k)
                   + self.H_shelf * (T_shelf_k - T_k))  # [XXX] units? - DRO

            # SOLID(IFYING) VIALS
            # heat capacity of solidifying vials
            # a vector of cps for all solidifying vials
            cp_sigma = (solid_fraction * cp_s
                        + (1-solid_fraction) * (sigma_k[solidMask] * (cp_i - cp_w) + cp_w))
            beta = depression * mass * cp_sigma

            deltaSigma = q_k[solidMask] / (alpha - beta / (1-sigma_k[solidMask])**2) * self.dt
            sigma_k[solidMask] = sigma_k[solidMask] + deltaSigma
            # assumption is that during solidification T = Teq
            T_k[solidMask] = T_eq - depression * (1/(1-sigma_k[solidMask]))

            # LIQUID VIALS
            # DRO: Leif, I've moved the dt out of the Hs because I associate Q with fluxes
            T_k[liquidMask] = T_k[liquidMask] + q_k[liquidMask] / hl * self.dt
            # a mask that is True where the temperature is below the equilibrium temperature
            superCooledMask = (T_k < T_eq)
            # a vial that is both liquid and supercooled can nucleate
            nucleationCandidatesMask = liquidMask & superCooledMask
            # the total number of nucleation candidates
            n_nucleationCandidates = np.sum(nucleationCandidatesMask)

            # nucleation probabilities for candidates
            P = np.zeros((self.N_vials_total,))
            diceRolls = np.zeros((self.N_vials_total,))

            P[nucleationCandidatesMask] = (kb * V * (T_eq - T_k[nucleationCandidatesMask])**b
                                           * self.dt)
            diceRolls[nucleationCandidatesMask] = np.random.rand(n_nucleationCandidates)

            # CONTROLLED NUCLEATION XXX

            # Nucleation
            nucleatedVialsMask = nucleationCandidatesMask & (diceRolls < P)

            stats['nucleationTimes'][nucleatedVialsMask] = t[k] + self.dt
            stats['nucleationTemperatures'][nucleatedVialsMask] = T_k[nucleatedVialsMask]

            q0 = (T_eq - T_k[nucleatedVialsMask]) * cp_solution * mass

            sigma_k[nucleatedVialsMask] = -q0 / (alpha - beta_solution)
            # assumption is that during solidification T = Teq
            T_k[nucleatedVialsMask] = T_eq - depression * (1/(1-sigma_k[nucleatedVialsMask]))

        self.stats = stats
        self._X = X  # store the state matrix
        self._t = t  # store the time vector

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
            t_nucleation = self.stats['nucleationTimes']

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
            T_nucleation = self.stats['nucleationTemperatures']

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
            t_solidification = self.stats['solidificationTimes']

        I_groups = self.getVialGroup(group)

        return t_solidification[I_groups]

    def solidCounter(self, threshold: float = 0.9) -> np.ndarray:
        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run before induction times can be extracted.")

        return np.sum(self.X_sigma > threshold, axis=0)

    def getVialGroup(self, group: Union[str, Sequence[str]] = 'all') -> np.ndarray:

        if isinstance(group, str):
            group = [group]

        _, VIAL_EXT = self._buildInteractionMatrices()

        myMask = np.zeros(self.N_vials_total, dtype=bool)
        for g in group:
            if g == 'corner':
                myMask = myMask | (VIAL_EXT == 2)
            elif g == 'edge':
                myMask = myMask | (VIAL_EXT == 1)
            elif g == 'center':
                myMask = myMask | (VIAL_EXT == 0)
            elif g == 'all':
                myMask = np.ones(self.N_vials_total, dtype=bool)
                break
            else:
                raise ValueError(f"Group {g} is not a known vial group.")

        return myMask

    def plot(self, what: str = 'temperature', kind: str = 'trajectories', group: Union[str, Sequence[str]] = 'all'):

        self.toDataframe()
        df = self.stats_df
        if group != 'all':
            df = df[df.group.isin(group)]

        if not kind.lower().startswith('traj'):
            df = df[df.variable.str.lower().str.contains(what.lower())]
            sns.catplot(data=df, hue='group', y='value', kind=kind, x='variable')

    def _sigmaCrossingIndices(self, threshold=0.9) -> Tuple[np.ndarray, np.ndarray]:

        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run before induction times can be extracted.")

        Indices = np.argmax(self.X_sigma > threshold, axis=1)
        # vials that never exceeded solidification threshold
        neverReached = ~np.any(self.X_sigma > threshold, axis=1)

        return Indices, neverReached

    def _buildInteractionMatrices(self) -> Tuple[csr_matrix, np.ndarray]:

        # This is a (n_x*n_y*n_z) x (n_x*n_y*n_z) matrix of interactions between vial pairs
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
        # pairs (1,2), (2,3), ..., (i,i+1) have horizontal interactions, except where i = n_x!
        # this means that in the IA matrix off-diagonal elements are 1
        # except where i%n_x==0 or j%n_x==0
        dx_pattern = np.ones((n_x*n_y-1,))
        # the vial at index position i+1 is at the edge -> one neighbor less
        idx_delete = [i for i in range(n_x*n_y-1) if (i+1) % n_x == 0]
        dx_pattern[idx_delete] = 0
        # we store this as a compressed sparse row (most efficient format?)
        DX = csr_matrix(np.diag(dx_pattern, k=1) + np.diag(dx_pattern, k=-1))

        # create interaction matrix for vertical (y-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,n_x+1), (2,n_x+2), ..., (i,n_x+i) have vertical interactions
        # for i in [1, n_x*(n_y-1)]
        # this means that in the IA matrix n_x-removed-off-diagonal elements are 1
        dy_pattern = np.ones((n_x*n_y - n_x,))
        DY = csr_matrix(np.diag(dy_pattern, k=n_x) + np.diag(dy_pattern, k=-n_x))

        # how many interactions does each vial have with other vials
        # diagflat because sum over a sparse matrix returns a np.matrix object, not an ndarray
        VIAL_INT = csr_matrix(np.diagflat(np.sum(DX + DY, axis=1)))

        # the overall interaction matrix is given by the sum of the interaction matrices DX and DY
        # minus the diagonal containing the sum of all vial-vial interactions
        interactionMatrix = DX + DY - VIAL_INT

        # at most, any cubic vial on a 2D shelf can have 4 interactions. 4 - VIAL_INT is the
        # number of external interactions (excl. the shelf)
        # is it worth storing this as a sparse matrix? - DRO XXX
        VIAL_EXT = 4*np.ones((n_x*n_y,)) - VIAL_INT.diagonal()

        return interactionMatrix, VIAL_EXT

    def _buildHeatflowMatrices(self):

        interactionMatrix, VIAL_EXT = self._buildInteractionMatrices()

        self.H_int = interactionMatrix * self.k['int'] * A

        self.H_ext = VIAL_EXT * self.k['ext'] * A

        if 's_sigma_rel' in self.k.keys():
            self.k['shelf'] = self.k['s0'] + np.random.normal(self.N_vials_total)
        else:
            self.k['shelf'] = self.k['s0']

        self.H_shelf = self.k['shelf'] * A  # either a scalar or a vector

    def toDataframe(self):

        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run before induction times can be extracted.")

        df = pd.DataFrame(self.stats)
        df.index.name = 'vial'
        df = df.reset_index()

        _, VIAL_EXT = self._buildInteractionMatrices()

        df['group'] = VIAL_EXT
        df.loc[df.group == 2, 'group'] = 'corner'
        df.loc[df.group == 1, 'group'] = 'edge'
        df.loc[df.group == 0, 'group'] = 'center'

        self.stats_df = df.melt(id_vars=['group', 'vial'])

    def __repr__(self) -> str:
        """ The string representation of the Snowfall class.

        Returns:
            str: The Snowfall class string representation giving some basic info.
        """

        return (f"Snowfall([N_vials: {self.N_vials}, Nrep: {self.Nrep}, "
                + f"dt: {self.dt}, pool_size: {self.poolsize}])")
