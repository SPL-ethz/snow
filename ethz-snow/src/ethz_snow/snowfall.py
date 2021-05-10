'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

from typing import List, Tuple
from operatingConditions import OperatingConditions
from scipy.sparse import csr_matrix

from constants import A

class Snowfall:
    def __init__(
        k: dict = {'int':20, 'ext': 20, 's0':20, 's_sigma_rel':0.1},
        N_vials: Tuple[int] = (7,7,1), # should this be part of operating conditions?
        Nrep: int = 5e3,
        dt: float = 2,
        pool_size: int = 12,
        opcond: OperatingConditions = OperatingConditions()
    ):

        if N_vials[2] > 1:
            raise NotImplementedError("Only 2D (shelf) models are implemented at this moment.")

    def run():
        
        # 1. build interaction matrices
        self._buildInteractionMatrices()
        
        # 2. fix k_shelf
        self._kShelf()

    def _buildInteractionMatrices():
        
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

        # at most, any cubic vial on a 2D shelf can have 4 interactions. 4 - VIAL_INT is the number of external interactions
        VIAL_EXT = csr_matrix(np.diag(4*np.ones((n_x*n_y,)))) - VIAL_INT

        # the overall interaction matrix is given by the sum of the interaction matrices DX and DY minus the diagonal containing the sum of all vial-vialinteractions
        interactionMatrix = DX + DY - VIAL_INT

        self.H_int = interactionMatrix * self.k['int'] * A * self.dt

        self.H_ext = VIAL_EXT * self.k['ext'] * A * self.dt

    def _kShelf():

        if 's_sigma_rel' in self.k.keys():

            self.k['shelf'] = self.k['s0'] + np.random.normal(np.prod(self.N_vials))

        else:
            
            self.k['shelf'] = self.k['s0']