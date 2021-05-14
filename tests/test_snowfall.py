'''
Created Date: Friday May 14th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

import pytest
from ethz_snow.snowfall import Snowfall
import numpy as np


def test_k_must_be_fine():
    with pytest.raises(TypeError):
        S = Snowfall(k=[1, 2, 3])

    # make sure k is checked w.r.t. its keys
    with pytest.raises(ValueError):
        S = Snowfall(k={'int': 1, 'ext': 1, 's0f': 1})


def test_opcond_must_be_operatingcondition():
    with pytest.raises(TypeError):
        S = Snowfall(opcond=dict(a=1))


def test_only2D_implemented():
    with pytest.raises(NotImplementedError):
        S = Snowfall(N_vials=(3, 3, 3))


@pytest.mark.parametrize('input', ['gibberish', '2random3', [1, 'asdf'], [0, 9], (-1, 2)])
def test_storeStatesMeaningless(input):
    with pytest.raises(ValueError):
        S = Snowfall(N_vials=(3, 3, 1), storeStates=input)


@pytest.mark.parametrize('input_result', [(None, 0), ('all', 9), ('corner', 4),
                                          ('edge', 4), ('center', 1), (['corner', 'edge'], 8),
                                          ('random4', 4), ('uniform+3', 3),
                                          (('center', 'corner_random_2'), 3)])
def test_storageMaskFunction(input_result):
    S = Snowfall(N_vials=(3, 3, 1), storeStates=input_result[0])

    assert np.sum(S._storageMask) == input_result[1]


@pytest.fixture(scope='module')
def S_331_all():
    S = Snowfall(N_vials=(3, 3, 1), storeStates='all')
    S.run()
    return S


@pytest.fixture(scope='module')
def S_331_edge():
    S = Snowfall(N_vials=(3, 3, 1), storeStates='edge')
    S.run()
    return S


@pytest.mark.parametrize('fun', ('nucleationTimes', 'nucleationTemperatures',
                         'solidificationTimes'))
def test_Statsconsistency(S_331_all, S_331_edge, fun):

    # the nucleation Temperatures can actually be a bit different, because the
    # fromStates method retrieves only the last Temperature before the nucleation
    # while the default method retrieves the stored, actual T of the vial
    assert abs(getattr(S_331_all, fun)().mean()
               - getattr(S_331_all, fun)(fromStates=True).mean()) < 0.1

    assert abs(getattr(S_331_edge, fun)(group='edge').mean()
               - getattr(S_331_edge, fun)(group='edge', fromStates=True).mean()) < 0.1
