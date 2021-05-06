'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

from typing import List
from operatingConditions import OperatingConditions
from scipy.sparse import csr_matrix

class Snowfall:
    def __init__(
        k: dict = {'int':20, 'ext': 20, 's0':20, 's_sigma_rel':0.1},
        N_vials: List[int] = [7,7,1], # should this be part of operating conditions?
        Nrep: int = 5e3,
        dt: float = 2,
        pool_size: int = 12,
        opcond: OperatingConditions = OperatingConditions()
    ):

    def run():
        # 1. build interaction matrices
        # 2. ...

    def _buildInteractionMatrices():
        
        # an interaction matrix is a (n_x*n_y*n_z) x (n_x*n_y*n_z) matrix of interactions between vial pairs

        n_x = self.N_vials[0] # the number of vials in the horizontal direction
        n_y = self.N_vials[1] # the number of vials in the vertical direction
        n_z = self.N_vials[2] # the number of vials in the upwards/downwards direction (perpendicular to the plane spanned by x and y)

        # create interaction matrix for horizontal (x-direction) interactions
        # matrix is 1 where two vials have a horizontal interaction and 0 otherwise
        # pairs (1,2), (2,3), ..., (i,i+1) have horizontal interactions, except where i = n_x!
        # this means that in the IA matrix off-diagonal elements are 1, except where i%n_x==0 or j%n_x==0
        dx_pattern = np.ones((n_x*n_y,))
        idx_delete = [i for i in range(n_x*n_y-1) if (i+1)%n_x == 0] # the vial at index position i+1 is at the edge -> one neighbor less
        dx_pattern[idx_delete] = 0
        DX = csr_matrix(np.diag(dx_pattern, k = 1) + np.diag(dx_pattern, k = -1)) # we store this as a compressed sparse row (most efficient format?)

        # create interaction matrix for horizontal (x-direction) interactions
        # MATLAB CODE
        # %y interaction (up and down)

        # dy = ones(x*y*z-x,1);

        # for u = 1:(x*y*z-x-1)
        #     if mod(u,(y*x)) == 0
        #         for k = 1:x
        #             dy(u+1-k) = 0;
        #         end
        #         u = u+x;
        #     end
        # end

        # %z interaction (front and back)

        # dz = ones(x*y*(z-1),1);

        # %Final matrix

        # INT = sparse(diag(d0) + diag(dx,1) + diag(dx,-1) + diag(dy,x) + diag(dy,-x));