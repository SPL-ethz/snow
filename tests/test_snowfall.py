"""Implement unit tests for Snowfall class."""
import pytest
from ethz_snow.snowfall import Snowfall
from ethz_snow.snowflake import Snowflake
from ethz_snow.operatingConditions import OperatingConditions

import numpy as np


def test_passSnowflakeArgs():
    """Test that Snowflake arguments are passed correctly."""
    S = Snowfall(N_vials=(3, 3, 1), dt=5, seed=121, solidificationThreshold=0.5)
    Sf = Snowflake()

    # seed input is deleted
    assert S.Sf_template.seed == Sf.seed
    assert S.Sf_template.dt == 5
    assert S.Sf_template.N_vials == (3, 3, 1)
    assert S.Sf_template.solidificationThreshold == 0.5


@pytest.fixture(scope="module")
def myOpcond():
    """Fixture to share operating conditions across tests."""
    cooling = {"start": 20, "end": -20, "rate": 1}
    myOpcond = OperatingConditions(t_tot=10000, cooling=cooling)

    return myOpcond


@pytest.fixture(scope="module")
def S_seq(myOpcond):
    """Fixture to only calculate S_seq once."""
    S_seq = Snowfall(N_vials=(3, 3, 1), dt=5, opcond=myOpcond, Nrep=5)
    S_seq.run(how="sequential")

    return S_seq


@pytest.mark.slow
def test_runTypesIdentical(S_seq, myOpcond):
    """Test that output is independent of parallelization."""
    S_async = Snowfall(N_vials=(3, 3, 1), dt=5, opcond=myOpcond, Nrep=5)
    S_async.run(how="async")
    t_nuc_async = np.concatenate(
        [S_async.stats[key]["t_nucleation"] for key in S_async.stats]
    )
    # need to sort because order of seeds might change in parallelization
    t_nuc_async.sort()

    S_sync = Snowfall(N_vials=(3, 3, 1), dt=5, opcond=myOpcond, Nrep=5)
    S_sync.run(how="sync")
    t_nuc_sync = np.concatenate(
        [S_sync.stats[key]["t_nucleation"] for key in S_sync.stats]
    )
    t_nuc_sync.sort()

    t_nuc_seq = np.concatenate(
        [S_seq.stats[key]["t_nucleation"] for key in S_seq.stats]
    )
    t_nuc_seq.sort()

    assert len(t_nuc_async) == 5 * 9
    assert all(t_nuc_async == t_nuc_sync)
    assert all(t_nuc_seq == t_nuc_sync)


@pytest.mark.slow
def test_statsComputation(S_seq):
    """Test that to_frame doesn't mix things up."""
    df = S_seq.to_frame()

    assert (
        df.loc[
            (df.seed == 0) & (df.vial == 1) & (df.variable == "T_nucleation"), "value"
        ].squeeze()
        == S_seq.stats[0]["T_nucleation"][1]
    )
