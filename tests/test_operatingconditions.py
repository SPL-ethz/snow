"""
Created Date: Wednesday May 19th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
"""
from _pytest.assertion import pytest_sessionfinish
import pytest
from ethz_snow.operatingConditions import OperatingConditions

import numpy as np


def test_holdingType():
    with pytest.raises(TypeError):
        OperatingConditions(holding=1)


@pytest.mark.parametrize(
    ("cooling", "holding"),
    [
        (dict(), dict()),
        ({"start": 10, "end": 1, "rate": 1}, dict()),
        ({"start": 10, "end": 1, "rate": 1}, {"duration": 10}),
        ({"start": 10, "end": 1}, {"duration": 10, "temp": 0}),
    ],
)
def test_requiredKeys(cooling, holding):
    with pytest.raises(ValueError):
        OperatingConditions(holding=holding, cooling=cooling)


@pytest.fixture(scope="module")
def myO():
    myO = OperatingConditions(
        t_tot=100,
        cooling=dict(start=20, end=-20, rate=1),
        holding=dict(duration=1, temp=0),
        controlledNucleation=False,
    )
    return myO


def test_cnt(myO):
    assert myO.cnt == np.inf
    myO.controlledNucleation = True
    assert myO.cnt == 21


def test_tempprofile(myO):
    T = myO.tempProfile(dt=1)

    assert T[0] == 20
    assert T[-1] == -20
    assert len(T) == 101
    assert np.sum(T == 0) == 2
    assert all(np.diff(T[:20]) == -1)
    assert all(T[41:] == -20)

    myO.holding["duration"] = 10
    T = myO.tempProfile(dt=3)
    assert all(T[7:11] == 0)
    assert T[11] != 0
    assert len(T) == 35
    assert all(np.diff(T[:7]) == -3)

    myO.cooling["rate"] = 2.5
    T = myO.tempProfile(dt=5)
    assert T[0] - T[1] == 12.5

    myO.cooling["rate"] = 1
    myO.holding = None
    myO.t_tot = 40
    T = myO.tempProfile(dt=1)
    assert len(T) == 41
    assert T[-1] == -20
