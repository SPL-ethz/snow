'''
Created Date: Friday May 14th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

import pytest
from ethz_snow.snowflake import Snowflake
import numpy as np


def test_k_must_be_fine():
    with pytest.raises(TypeError):
        S = Snowflake(k=[1, 2, 3])

    # make sure k is checked w.r.t. its keys
    with pytest.raises(ValueError):
        S = Snowflake(k={'int': 1, 'ext': 1, 's0f': 1})


def test_opcond_must_be_operatingcondition():
    with pytest.raises(TypeError):
        S = Snowflake(opcond=dict(a=1))


def test_only2D_implemented():
    with pytest.raises(NotImplementedError):
        S = Snowflake(N_vials=(3, 3, 3))


@pytest.mark.parametrize('input', ['gibberish', '2random3', [1, 'asdf'], [0, 9], (-1, 2)])
def test_storeStatesMeaningless(input):
    with pytest.raises(ValueError):
        S = Snowflake(N_vials=(3, 3, 1), storeStates=input)


@pytest.mark.parametrize('input_result', [(None, 0), ('all', 9), ('corner', 4),
                                          ('edge', 4), ('center', 1), (['corner', 'edge'], 8),
                                          ('random4', 4), ('uniform+3', 3),
                                          (('center', 'corner_random_2'), 3)])
def test_storageMaskFunction(input_result):
    S = Snowflake(N_vials=(3, 3, 1), storeStates=input_result[0])

    assert np.sum(S._storageMask) == input_result[1]


@pytest.fixture(scope='module')
def fakeS():
    S = Snowflake(N_vials=(2, 2, 1), storeStates='all')
    S._t = np.array([0, 1, 2, 3])
    S._X = np.concatenate([np.zeros((4, 4)),
                           np.array([[0, 0, 0.5, 0.95], [0, 0.7, 0.75, 1],
                                    [0, 0, 0, 0], [0, 1, 1, 1]])], axis=0)
    S.stats['t_solidification'] = np.array([3, 3, np.nan, 1])
    S.stats['t_nucleation'] = np.array([2, 1, np.nan, 1])
    S.stats['T_nucleation'] = np.array([0, 0, np.nan, 0])

    return S


def test_sigmaCounter(fakeS):
    assert fakeS.sigmaCounter(time=1) == 1
    assert all(fakeS.sigmaCounter(time=[0, 1, 2, 3]) == np.array([0, 1, 1, 3]))
    assert fakeS.sigmaCounter(time=2) == fakeS.sigmaCounter(time=2, fromStates=True)
    assert fakeS.sigmaCounter(time=2, threshold=0.4, fromStates=True) == 3


def test_getvialgroup(S_331_all):

    assert np.where(S_331_all.getVialGroup('center'))[0] == 4
    assert all(np.where(S_331_all.getVialGroup(['center', 'corner']))[0] == [0, 2, 4, 6, 8])

    with pytest.raises(ValueError):
        S_331_all.getVialGroup('rubbish')


def test_toDataframe(fakeS):

    stats_df, traj_df = fakeS.toDataframe(n_timeSteps=4)

    assert all([col in stats_df.columns for col in ['group', 'vial', 'variable', 'value']])
    assert stats_df.shape[0] == 12
    assert all(stats_df.group == 'corner')

    assert traj_df.loc[(traj_df.vial == 0) & (traj_df.state == 'sigma') & (traj_df.Time == 3), 'value'].item() == 0.95


@pytest.fixture(scope='module')
def S_331_all():
    S = Snowflake(N_vials=(3, 3, 1), storeStates='all')
    S.run()
    return S


@pytest.fixture(scope='module')
def S_331_edge():
    S = Snowflake(N_vials=(3, 3, 1), storeStates='edge')
    S.run()
    return S

@pytest.mark.slow
@pytest.mark.parametrize('fun', ('nucleationTimes', 'nucleationTemperatures',
                         'solidificationTimes'))
def test_statsConsistency(S_331_all, S_331_edge, fun):

    # the nucleation Temperatures can actually be a bit different, because the
    # fromStates method retrieves only the last Temperature before the nucleation
    # while the default method retrieves the stored, actual T of the vial
    assert abs(getattr(S_331_all, fun)().mean()
               - getattr(S_331_all, fun)(fromStates=True).mean()) < 0.1

    assert abs(getattr(S_331_edge, fun)(group='edge').mean()
               - getattr(S_331_edge, fun)(group='edge', fromStates=True).mean()) < 0.1
