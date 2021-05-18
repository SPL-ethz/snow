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

import time

HEATFLOW_REQUIREDKEYS = ('int', 'ext', 's0')
VIAL_GROUPS = ('corner', 'edge', 'center', 'all')


class Snowflake:
    def __init__(
        self,
        k: dict = {'int': 20, 'ext': 20, 's0': 20, 's_sigma_rel': 0.1},
        N_vials: Tuple[int, int, int] = (7, 7, 1),  # should this be part of operating conditions?
        storeStates: Optional[Union[str, Sequence[str], Sequence[int]]] = None,
        solidificationThreshold: float = 0.9,
        dt: float = 2,
        seed: int = 2021,
        opcond: OperatingConditions = OperatingConditions()
    ):

        if not isinstance(k, dict):
            raise TypeError(f"Input k must be of type dict. Was given {type(k)}.")
        elif not all([key in k.keys() for key in HEATFLOW_REQUIREDKEYS]):
            raise ValueError((f"A required key was missing from dictionary k, specifically "
                             + f"{set(HEATFLOW_REQUIREDKEYS) - set(k.keys())}."))
        else:
            self.k = k

        if isinstance(N_vials, list):
            # we exepct this to be immutable
            N_vials = tuple(N_vials)

        if N_vials[2] > 1:
            raise NotImplementedError("Only 2D (shelf) models are implemented at this moment.")
        else:
            self.N_vials = N_vials

        self.dt = dt

        # store the seed to look it up if need be
        self.seed = seed

        if not isinstance(opcond, OperatingConditions):
            raise TypeError("Input opcond must be an instance of class OperatingConditions.")
        else:
            self.opcond = opcond

        if isinstance(storeStates, (list, tuple)):
            if all([isinstance(x, int) for x in storeStates]):
                if (any(np.array(storeStates) > self.N_vials_total - 1)
                        or any(np.array(storeStates) < 0)):
                    raise ValueError(f"Entries in storeStates must be >0 "
                                     + f"and <{self.N_vials_total - 1}.")
                storageMask = np.zeros(self.N_vials_total, dtype=bool)
                storageMask[storeStates] = True
            elif all([isinstance(x, str) for x in storeStates]):
                storageMasks = list(map(lambda x: self._interpretStorageString(x), storeStates))
                storageMask = np.logical_or.reduce(storageMasks)
            else:
                raise ValueError("storeStates must be a sequence of all int or all str.")
        elif isinstance(storeStates, str):
            storageMask = self._interpretStorageString(storeStates)
        elif storeStates is None:
            storageMask = np.zeros(self.N_vials_total, dtype=bool)

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
    def simulationStatus(self):
        if (self._simulationStatus == 1) or (self._X is not None):
            self._simulationStatus = 1

        return self._simulationStatus

    @property
    def H_int(self):
        if (self._H_int is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
        return self._H_int

    @property
    def H_ext(self):
        if (self._H_ext is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
        return self._H_ext

    @property
    def H_shelf(self):
        if (self._H_shelf is None) or self.N_vials != self._NvialsUsed:
            self._buildHeatflowMatrices()
        return self._H_shelf

    @property
    def N_vials_total(self):
        return int(np.prod(self.N_vials))

    @property
    def X_T(self):
        return self._X[:int(np.sum(self._storageMask)), :]

    @property
    def X_sigma(self):
        return self._X[int(np.sum(self._storageMask)):, :]

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, value):
        # best practice is now to store the rng in an object and pass it around
        # https://towardsdatascience.com/stop-using-numpy-random-seed-581a9972805f
        self._rng = np.random.default_rng(value)
        # store the seed to look it up if need be
        self._seed = value

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
                I_toStore = self._rng.choice(I_candidates, size=howMany, replace=False)
            elif 'uniform' in myString:
                stepSize = int(np.ceil(len(I_candidates)/howMany))
                I_fromCandidates = np.arange(0, len(I_candidates), stepSize, dtype=int)
                I_toStore = I_candidates[I_fromCandidates]
            storageMask = np.zeros(self.N_vials_total, dtype=bool)
            storageMask[I_toStore] = True

        return storageMask

    def run(self):
        # tic = time.perf_counter()
        
        # clean up any potential old simulations
        self._X = None
        self._t = None

        N_timeSteps = int(np.ceil(self.opcond.t_tot / self.dt))+1  # the total number of timesteps
        # toc1 = time.perf_counter()
        # print(toc1-tic)
        # toc2 = time.perf_counter()
        # print(f"BuildHeatflow: {toc2-toc1}")
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
            t_nucleation=np.full(self.N_vials_total, np.nan),
            T_nucleation=np.full(self.N_vials_total, np.nan),
            t_solidification=np.full(self.N_vials_total, np.nan))

        k_CN = np.argmax(t >= self.opcond.cnt)
        # toc3 = time.perf_counter()
        # print(f"Tempprofile: {toc3-toc2}")
        # 5. Iterate over time steps
        for k in np.arange(N_timeSteps):

            T_shelf_k = T_shelf[k]
            T_ext_k = T_ext[k]
            # toc1_l = time.perf_counter()
            x_k = np.concatenate([T_k, sigma_k])
            X[:, k] = x_k[np.tile(self._storageMask, 2)]

            # a mask that is True where a vial is still fully liquid
            liquidMask = (sigma_k == 0)
            solidMask = ~liquidMask  # for convenience
            solidifiedMask = ((sigma_k > self.solidificationThreshold)
                              & np.isnan(stats['t_solidification']))
            stats['t_solidification'][solidifiedMask] = (t[k]
                                                         - stats['t_nucleation'][solidifiedMask])

            # toc2_l = time.perf_counter()
            # print(f"Masks: {toc2_l-toc1_l:4.2e}")
            # calculate heatflows for all vials
            q_k = (self.H_int @ T_k
                   + self.H_ext * (T_ext_k - T_k)
                   + self.H_shelf * (T_shelf_k - T_k))  # [XXX] units? - DRO
            # toc2b_l = time.perf_counter()
            # print(f"q: {toc2b_l-toc2_l:4.2e}")
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
            # toc3_l = time.perf_counter()
            # print(f"Time steps: {toc3_l-toc2b_l:4.2e}")

            # nucleation probabilities for candidates
            P = np.zeros((self.N_vials_total,))
            diceRolls = np.zeros((self.N_vials_total,))

            P[nucleationCandidatesMask] = (kb * V * (T_eq - T_k[nucleationCandidatesMask])**b
                                           * self.dt)
            # toc4_l = time.perf_counter()
            # print(f"Probabilities: {toc4_l-toc3_l:4.2e}")
            # when we reach the timepoint of controlled nucleation
            # all vials (that thermodynamically can) nucleate
            if k == k_CN:
                P.fill(1)
            diceRolls[nucleationCandidatesMask] = self._rng.random(n_nucleationCandidates)

            # CONTROLLED NUCLEATION XXX

            # Nucleation
            nucleatedVialsMask = nucleationCandidatesMask & (diceRolls < P)

            stats['t_nucleation'][nucleatedVialsMask] = t[k] + self.dt
            stats['T_nucleation'][nucleatedVialsMask] = T_k[nucleatedVialsMask]
            # toc5_l = time.perf_counter()
            # print(f"Dice rolls: {toc5_l-toc4_l:4.2e}")

            q0 = (T_eq - T_k[nucleatedVialsMask]) * cp_solution * mass

            sigma_k[nucleatedVialsMask] = -q0 / (alpha - beta_solution)
            # assumption is that during solidification T = Teq
            T_k[nucleatedVialsMask] = T_eq - depression * (1/(1-sigma_k[nucleatedVialsMask]))
            # sys.exit()

        # toc4 = time.perf_counter()
        # print(f"For loop: {toc4-toc3:4.2e}")
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

        self._NvialsUsed = self.N_vials

        interactionMatrix, VIAL_EXT = self._buildInteractionMatrices()

        self._H_int = interactionMatrix * self.k['int'] * A

        self._H_ext = VIAL_EXT * self.k['ext'] * A

        if 's_sigma_rel' in self.k.keys():
            self.k['shelf'] = self.k['s0'] + self._rng.normal(size=self.N_vials_total)
        else:
            self.k['shelf'] = self.k['s0']

        self._H_shelf = self.k['shelf'] * A  # either a scalar or a vector

    def toDataframe(self, n_timeSteps=250) -> Tuple[np.ndarray, Optional[np.ndarray]]:

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

        stats_df = df.melt(id_vars=['group', 'vial'])

        n_storedStates = np.sum(self._storageMask)
        if n_storedStates > 0:
            # to reduce computational cost we limit number of timepoints to n_timeSteps
            df = pd.DataFrame(self._X[:, ::int(self._X.shape[1]/(n_timeSteps - 1))])
            df.columns = self._t[::int(self._X.shape[1]/(n_timeSteps - 1))]
            df.columns.name = 'Time'

            df['state'] = np.concatenate([np.repeat('temperature', n_storedStates),
                                          np.repeat('sigma', n_storedStates)])
            df['vial'] = np.tile(np.where(self._storageMask)[0], 2)
            df['group'] = np.tile(VIAL_EXT[self._storageMask], 2)

            df.loc[df.group == 2, 'group'] = 'corner'
            df.loc[df.group == 1, 'group'] = 'edge'
            df.loc[df.group == 0, 'group'] = 'center'

            traj_df = df.melt(id_vars=['group', 'vial', 'state'])
        else:
            traj_df = None

        return stats_df, traj_df

    def __repr__(self) -> str:
        """ The string representation of the Snowflake class.

        Returns:
            str: The Snowflake class string representation giving some basic info.
        """

        return (f"Snowflake([N_vials: {self.N_vials}, "
                + f"dt: {self.dt}, seed: {self.seed}])")
