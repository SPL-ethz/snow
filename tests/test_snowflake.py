"""Implement unit tests for Snowflake class."""
import pytest
from ethz_snow.snowflake import Snowflake
import numpy as np


def test_k_must_be_fine():
    """Test that k is checked properly."""
    with pytest.raises(TypeError):
        S = Snowflake(k=[1, 2, 3])

    # make sure k is checked w.r.t. its keys
    with pytest.raises(ValueError):
        S = Snowflake(k={"int": 1, "ext": 1, "s0f": 1})


def test_Constant_kShelf_fine():
    """Check that constant kshelf is respected."""
    d = {"int": 20, "ext": 20, "s0": 20, "s_sigma_rel": 0}
    S = Snowflake(storeStates="all", N_vials=(7, 7, 1), k=d)
    S._buildHeatflowMatrices()

    k_shelf = S.k["shelf"]

    assert isinstance(k_shelf, (int, float))

    assert k_shelf == 20


def test_opcond_must_be_operatingcondition():
    """Check that opcond type is verified."""
    with pytest.raises(TypeError):
        S = Snowflake(opcond=dict(a=1))


def test_only2D_implemented():
    """Ensure no 3D model is being asked."""
    with pytest.raises(NotImplementedError):
        S = Snowflake(N_vials=(3, 3, 3))


@pytest.mark.parametrize(
    "input", ["gibberish", "2random3", [1, "asdf"], [0, 9], (-1, 2)]
)
def test_storeStatesMeaningless(input):
    """Ensure invalid storeStates are rejected."""
    with pytest.raises(ValueError):
        S = Snowflake(N_vials=(3, 3, 1), storeStates=input)


@pytest.mark.parametrize(
    "input_result",
    [
        (None, 0),
        ("all", 9),
        ("corner", 4),
        ("edge", 4),
        ("center", 1),
        (["corner", "edge"], 8),
        ("random4", 4),
        ("uniform+3", 3),
        (("center", "corner_random_2"), 3),
    ],
)
def test_storageMaskFunction(input_result):
    """Ensure valid storeStates are interpreted correctly."""
    S = Snowflake(N_vials=(3, 3, 1), storeStates=input_result[0])

    assert np.sum(S._storageMask) == input_result[1]


@pytest.fixture(scope="module")
def fakeS():
    """Fixture to share a fake Snowflake across tests.
    
    Fake in that states are manually provided to avoid rng issues.
    """
    S = Snowflake(N_vials=(2, 2, 1), storeStates="all")
    S._t = np.array([0, 1, 2, 3])
    S._X = np.concatenate(
        [
            np.zeros((4, 4)),
            np.array(
                [[0, 0, 0.5, 0.95], [0, 0.7, 0.75, 1], [0, 0, 0, 0], [0, 1, 1, 1]]
            ),
        ],
        axis=0,
    )
    S.stats["t_solidification"] = np.array([3, 3, np.nan, 1])
    S.stats["t_nucleation"] = np.array([2, 1, np.nan, 1])
    S.stats["T_nucleation"] = np.array([0, 0, np.nan, 0])

    return S


def test_sigmaCounter(fakeS):
    """Test that sigmas are calculated correctly."""
    assert fakeS.sigmaCounter(time=1) == 1
    assert all(fakeS.sigmaCounter(time=[0, 1, 2, 3]) == np.array([0, 1, 1, 3]))
    assert fakeS.sigmaCounter(time=2) == fakeS.sigmaCounter(time=2, fromStates=True)
    assert fakeS.sigmaCounter(time=2, threshold=0.4, fromStates=True) == 3


def test_getvialgroup(S_331_all):
    """Test that vialgrouping masks are correctly computed."""
    assert np.where(S_331_all.getVialGroup("center"))[0] == 4
    assert all(
        np.where(S_331_all.getVialGroup(["center", "corner"]))[0] == [0, 2, 4, 6, 8]
    )

    with pytest.raises(ValueError):
        S_331_all.getVialGroup("rubbish")


def test_to_frame(fakeS):
    """Test that conversion to dataframe works as intended."""
    stats_df, traj_df = fakeS.to_frame(n_timeSteps=4)

    assert all(
        [col in stats_df.columns for col in ["group", "vial", "variable", "value"]]
    )
    assert stats_df.shape[0] == 12
    assert all(stats_df.group == "corner")

    assert (
        traj_df.loc[
            (traj_df.vial == 0) & (traj_df.state == "sigma") & (traj_df.Time == 3),
            "value",
        ].item()
        == 0.95
    )


@pytest.fixture(scope="module")
def S_331_all():
    """Fixture to share Snowflake across tests."""
    S = Snowflake(N_vials=(3, 3, 1), storeStates="all")
    S.run()
    return S


@pytest.fixture(scope="module")
def S_331_edge():
    """Fixture to share Snowflake across tests."""
    S = Snowflake(N_vials=(3, 3, 1), storeStates="edge")
    S.run()
    return S


@pytest.mark.slow
@pytest.mark.parametrize(
    "fun", ("nucleationTimes", "nucleationTemperatures", "solidificationTimes")
)
def test_statsConsistency(S_331_all, S_331_edge, fun):
    """Ensure statistics are calculated correctly.

    It should not matter whether we are using fromStates or not and whether we are
    storing all the states or not. Note that the nucleation Temperatures can 
    actually vary a bit, because the fromStates method retrieves only the last 
    temperature before the nucleation while the default method retrieves the stored, 
    _actual_ T of the vial at the time of nucleation.
    """
    #
    assert (
        abs(
            getattr(S_331_all, fun)().mean()
            - getattr(S_331_all, fun)(fromStates=True).mean()
        )
        < 0.1
    )

    assert (
        abs(
            getattr(S_331_edge, fun)(group="edge").mean()
            - getattr(S_331_edge, fun)(group="edge", fromStates=True).mean()
        )
        < 0.1
    )


def test_interactionMatrix():
    """Ensure interaction matrices are computed correctly."""
    S = Snowflake(N_vials=(3, 3, 1))
    mint, mext = S._buildInteractionMatrices()

    assert mint.shape == (9, 9)
    assert mint.nnz == 33
    assert mint[3, 0] == 1
    assert mint[0, 3] == 1
    assert mint[8, 5] == 1
    assert mint[4, 4] == -4
    assert (np.sum(mint, axis=0) == 0).all()
    assert (np.sum(mint, axis=1) == 0).all()

    assert (4 + mint.diagonal() == mext).all()
