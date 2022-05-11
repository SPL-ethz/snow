"""Implement unit tests for operatingConditions class."""
from _pytest.assertion import pytest_sessionfinish
import pytest
from ethz_snow.operatingConditions import OperatingConditions

import numpy as np
import warnings


def test_holdingType():
    """Make sure holding type is checked."""
    with pytest.raises(TypeError):
        OperatingConditions(holding=1)


@pytest.mark.parametrize(
    ("cooling", "holding"),
    [
        (dict(), dict()),
        ({"start": 10, "end": 1, "rate": 1}, dict()),
        ({"start": 10, "end": 1, "rate": 1}, {"duration": 10}),
        ({"start": 10, "end": 1}, {"duration": 10, "temp": 0}),
        ({"start": 10, "end": 1}, [{"duration": 10, "temp": 0}, {"duration": 3}]),
        ({"start": 20, "end": 10, "rate": 0}, None),
        ({"start": 20, "end": 20, "rate": 0}, {"duration": 10, "temp": 0}),
    ],
)
def test_requiredKeys(cooling, holding):
    """Make sure inputs are tested for required keys."""
    with pytest.raises(ValueError):
        OperatingConditions(holding=holding, cooling=cooling)


@pytest.mark.parametrize(
    ("cooling", "holding"),
    [
        ({"start": 10, "end": 1, "rate": 1}, None),
        ({"rate": 1e-16, "start": -20, "end": -20}, {"duration": 10, "temp": 0}),
    ],
)
def test_impliedTimesvsTotal(cooling, holding):
    with pytest.warns(UserWarning):
        OperatingConditions(
            t_tot=1,
            cooling=cooling,
            holding=holding,
        )


@pytest.fixture(scope="module")
def myO():
    """Fixture to share operatingconditions across tests."""
    myO = OperatingConditions(
        t_tot=100,
        cooling=dict(start=20, end=-20, rate=1),
        holding=dict(duration=1, temp=0),
        cnTemp=None,
    )
    return myO


@pytest.fixture(scope="module")
def myO_mh():
    """Fixture to share multihold operatingconditions across tests."""
    myO_mh = OperatingConditions(
        t_tot=500,
        cooling=dict(start=20, end=-20, rate=1),
        holding=[dict(duration=69, temp=-10), dict(duration=1, temp=0)],
        cnTemp=-10,
    )
    # cnt should be at 30 (time to cool down) +1 (first hold)+69 (second hold)=100
    return myO_mh


def test_cnt(myO, myO_mh):
    """Test computation of controlled nuc. time."""
    # if cn is false, cnt is inf
    assert myO.cnt == np.inf

    myO.cnTemp = myO.holding[0]["temp"]
    assert myO.cnt == 21

    assert myO_mh.cnt == 100


def test_tempprofile(myO, myO_mh):
    """Test that tempprofiles are calculated correctly."""
    T = myO.tempProfile(dt=1)

    assert T[0] == 20
    assert T[-1] == -20
    assert len(T) == 101
    assert np.sum(T == 0) == 2
    assert all(np.diff(T[:20]) == -1)
    assert all(T[41:] == -20)

    myO.holding[0]["duration"] = 10
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

    T_mh = myO_mh.tempProfile(1)

    assert len(T_mh) == 501
    assert all(T_mh[20:22] == 0)
    assert all(T_mh[31:101] == -10)
