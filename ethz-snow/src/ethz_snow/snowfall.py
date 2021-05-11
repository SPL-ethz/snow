'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

from operatingConditions import OperatingConditions
from scipy.sparse import csr_matrix

from constants import A, hl, T_eq, kb, V, b, cp_solution, mass, solid_fraction, cp_s, depression

from typing import List, Tuple, Union, Sequence

HEATFLOW_REQUIREDKEYS = ('int','ext','s0')
class Snowfall:
    def __init__(
        self,
        k: dict = {'int':20, 'ext': 20, 's0':20, 's_sigma_rel':0.1},
        N_vials: Tuple[int] = (7,7,1), # should this be part of operating conditions?
        Nrep: int = 5e3,
        dt: float = 2,
        pool_size: int = 12,
        opcond: OperatingConditions = OperatingConditions()
    ):

        if not isinstance(k, dict):
            raise TypeError(f"Input k must be of type dict. Was given {type(k)}.")
        elif not all([key in k.keys() for key in HEATFLOW_REQUIREDKEYS]):
            raise ValueError(f"A required key was missing from dictionary k, specifically {set(HEATFLOW_REQUIREDKEYS) - set(k.keys())}.")
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
        return np.prod(self.N_vials)

    @property
    def X_T(self):
        return self._X[:self.N_vials_total,:]

    @property
    def X_sigma(self):
        return self._X[self.N_vials_total:,:]

    def run():
        # clean up any potential old simulations
        self._X = None
        self._t = None

        N_timeSteps = np.ceil(self.opcond.t_tot / self.dt)+1 # the total number of timesteps

        # 1. build interaction matrices
        self._buildHeatflowMatrices()

        # 2. Obtain external temperature profile over time
        T_shelf = self.opcond.tempProfile(self.dt)
        T_ext = T_shelf # need to make a switch so that this can be decoupled - DRO XX

        # 3. Initial state of the system
        T_k = np.ones((self.N_vials_total,)) * self.opcond.cooling['start'] # [C] initial temperature
        sigma_k = np.zeros((self.N_vials_total,)) # state of the vials (0 = liq, 1 = completely frozen)
        T_shelf_k = T_shelf[0]
        T_ext_k = T_text[0]

        # 4. Pre-allocate memory
        X = np.zeros((2*self.N_vials_total, N_timeSteps)) # our state matrix containing the entire history of the system
        t = np.arange(0,N_timeSteps * self.dt, self.dt)

        # 5. Iterate over time steps
        for k in np.arange(N_timeSteps):

            X[:,k] = np.concatenate([T_k, sigma_k])

            liquidMask = (sigma_k == 0) # a mask that is True where a vial is still fully liquid
            solidMask = ~liquidMask # for convenience
            superCooledMask = (T_k < T_eq) # a mask that is True where the temperature is below the equilibrium temperature
            nucleationCandidatesMask = liquidMask & superCooledMask
            n_nucleationCandidates = np.sum(nucleationCandidatesMask)

            # calculate heatflows for all vials
            q_k = self.H_int @ T_k + self.H_ext * (T_ext_k - T_k) + self.H_shelf * (T_shelf_k - T_k) # [XXX] units? - DRO

            ## SOLID(IFYING) VIALS
            # heat capacity of solidifying vials
            cp_sigma = solid_fraction * cp_s + (1-solid_fraction) * (sigma_k[solidMask] * (cp_i - cp_w) + cp_w)

            ## LIQUID VIALS
            T_k[liquidMask] = T_k[liquidMask] + q_k[liquidMask] / hl

            # nucleation probabilities for candidates
            P = kb * V *(T_eq - T_k(nucleationCandidatesMask))**b * self.dt
            diceRolls = np.random.rand(n_nucleationCandidates)

            # CONTROLLED NUCLEATION XXX

            # Nucleation
            nucleatedVialsMask = nucleationCandidatesMask & (diceRolls < P)
            q0 = (T_eq - T_k[nucleatedVialsMask]) * cp_solution * mass 

        self._X = X # store the state matrix
        self._t = t

    def nucleationTimes(self):
        
        I_nucleation, lateBloomers = self._sigmaCrossingIndices(threshold = 0)

        # nucleation times for all nucleated vials
        t_nucleation = self._t[I_nucleation].astype(float) # need to make sure this is float so no problems arise later

        # non-nucleated vials are set to NaN
        t_nucleation[lateBloomers] = np.nan

        return t_nucleation

    def nucleationTemperatures(self) -> np.ndarray:
        
        I_nucleation, lateBloomers = self._sigmaCrossingIndices(threshold = 0)

        T_nucleation = np.array([self.X_T[i,I] for i,I in zip(range(self.N_vials_total), I_nucleation)]).astype(float) # this should always be float, but just to be sure

        # non-nucleated vials are set to NaN
        T_nucleation[lateBloomers] = np.nan

        return T_nucleation

    def solidificationTimes(self) -> np.ndarray:

        t_nucleation = self.nucleationTimes()

        I_solidification, neverGrownUp = self._sigmaCrossingIndices(threshold = 0.9)

        # nucleation times for all nucleated vials
        t_solidified = self._t[I_solidification].astype(float) # need to make sure this is float so no problems arise later

        # never fully solidified vials are set to NaN
        t_solidified[neverGrownUp] = np.nan

        # solidification is the difference between solidified and nucleated times
        t_solidification = t_solidified - t_nucleation

        return t_solidification

    def getVialGroup(self, group: Union[str, Sequence[str]]) -> np.ndarray:
        
        if isinstance(group, str):
            group = [group]

        _ , VIAL_EXT = self._buildInteractionMatrices()
        
        myMask = np.zeros(self.N_vials_total, dtype = bool)
        for g in group:
            if g == 'corner':
                myMask = myMask | (VIAL_EXT == 2)
            elif g == 'edge':
                myMask = myMask | (VIAL_EXT == 1)
            elif g == 'center':
                myMask = myMask | (VIAL_EXT == 0)

        return myMask

    def _sigmaCrossingIndices(self, thresold = 0.9) -> Tuple[np.ndarray, np.ndarray]:

        if self.simulationStatus == 0:
            raise ValueError("Simulation needs to be run before induction times can be extracted.")

        I= np.argmax(self.X_sigma > threshold, axis = 1)
        neverReached = ~np.any(self.X_sigma > threshold, axis = 1) # vials that never exceeded solidification threshold

        return I, neverReached

    def _buildInteractionMatrices(self):
        
        # an interaction matrix is a (n_x*n_y*n_z) x (n_x*n_y*n_z) matrix of interactions between vial pairs
        # Please note that, on any given shelf, we index vials this way:
        # [[1, 2, 3, 4, ..., n_x],
        # [n_x+1, n_x+2, ..., 2*n_x],
        # [...],
        # [(n_y-1)*n_x+1, ... , n_y*n_x]]
        # If one thinks of a shelf as an array, we hence use 'row-major' ordering. This is the same as numpy (most of the time), but different from Matlab.
        # These matrices can then be used as a sort of mask to overlay on, e.g., temperatures to calculate the correct driving forces for each vial/pair.

        n_x = self.N_vials[0] # the number of vials in the horizontal direction
        n_y = self.N_vials[1] # the number of vials in the vertical direction
        n_z = self.N_vials[2] # the number of vials in the upwards/downwards direction (perpendicular to the plane spanned by x and y)

        # create interaction matrix for horizontal (x-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,2), (2,3), ..., (i,i+1) have horizontal interactions, except where i = n_x!
        # this means that in the IA matrix off-diagonal elements are 1, except where i%n_x==0 or j%n_x==0
        dx_pattern = np.ones((n_x*n_y-1,))
        idx_delete = [i for i in range(n_x*n_y-1) if (i+1)%n_x == 0] # the vial at index position i+1 is at the edge -> one neighbor less
        dx_pattern[idx_delete] = 0
        DX = csr_matrix(np.diag(dx_pattern, k = 1) + np.diag(dx_pattern, k = -1)) # we store this as a compressed sparse row (most efficient format?)

        # create interaction matrix for vertical (y-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,n_x+1), (2,n_x+2), ..., (i,n_x+i) have vertical interactions for i in [1, n_x*(n_y-1)]
        # this means that in the IA matrix n_x-removed-off-diagonal elements are 1
        dy_pattern = np.ones((n_x*n_y - n_x,))
        DY = csr_matrix(np.diag(dy_pattern, k = n_x) + np.diag(dy_pattern, k = -n_x))

        # how many interactions does each vial have with other vials
        VIAL_INT = csr_matrix(np.diagflat(np.sum(DX + DY, axis = 1))) # diagflat because sum over a sparse matrix returns a np.matrix object, not an ndarray

        # the overall interaction matrix is given by the sum of the interaction matrices DX and DY minus the diagonal containing the sum of all vial-vial interactions
        interactionMatrix = DX + DY - VIAL_INT

        # at most, any cubic vial on a 2D shelf can have 4 interactions. 4 - VIAL_INT is the number of external interactions (excl. the shelf)
        VIAL_EXT = 4*np.ones((n_x*n_y,)) - VIAL_INT.diagonal() # is it worth storing this as a sparse matrix? - DRO XXX

        return interactionMatrix, VIAL_EXT

    def _buildHeatflowMatrices(self):

        interactionMatrix, VIAL_EXT = self._buildInteractionMatrices()

        self.H_int = interactionMatrix * self.k['int'] * A * self.dt

        self.H_ext = VIAL_EXT * self.k['ext'] * A * self.dt

        self._kShelf()
        self.H_shelf = self.k['shelf'] * A * self.dt # either a scalar or a vector

    def _kShelf(self):

        if 's_sigma_rel' in self.k.keys():

            self.k['shelf'] = self.k['s0'] + np.random.normal(np.prod(self.N_vials))

        else:
            
            self.k['shelf'] = self.k['s0']