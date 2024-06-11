"""Implement unit tests for Snowflake class."""

import os
import pytest
from ethz_snow.operatingConditions import OperatingConditions
from ethz_snow.snowflake import Snowflake
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.parametrize("inits", [[1, 2, 3], "something"])
def test_initialStates_fine(inits):
    with pytest.raises(ValueError):
        S = Snowflake(initialStates=inits)


@pytest.mark.parametrize(
    "inits", [dict(sigma=0.5, temp=20), dict(sigma=None, temp=(1, 2, 3))]
)
def test_sigmaInit_notImplemented(inits):
    with pytest.raises(NotImplementedError):
        S = Snowflake(initialStates=inits)


def test_k_must_be_fine():
    """Test that k is checked properly."""
    with pytest.raises(TypeError):
        S = Snowflake(k=[1, 2, 3])

    # make sure k is checked w.r.t. its keys
    with pytest.raises(ValueError):
        S = Snowflake(k={"int": 1, "ext": 1, "s0f": 1})

    # in case of 3D s0 is not required
    S = Snowflake(N_vials=(3, 3, 3), k=dict(int=1, ext=1))


@pytest.mark.parametrize("initIce", ["nondirect", dict(initIce=1), 3])
def test_initIce_fine(initIce):
    with pytest.raises(ValueError):
        S = Snowflake(initIce=initIce)


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


@pytest.mark.parametrize(
    "input",
    [
        (1, 1),
        (1,),
        (3, 3, 3, 3),
    ],
)
def test_N_vials_dims(input):
    """N_vials must be a tuple with length 3."""
    with pytest.raises(ValueError):
        S = Snowflake(N_vials=input)


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
        ("core", 1),
        (["corner", "edge"], 8),
        ("random4", 4),
        ("uniform+3", 3),
        (("core", "corner_random_2"), 3),
    ],
)
def test_storageMaskFunction_square(input_result):
    """Ensure valid storeStates are interpreted correctly."""
    config = THIS_DIR + "/data/square.yaml"
    S = Snowflake(N_vials=(3, 3, 1), storeStates=input_result[0], configPath=config)

    assert np.sum(S._storageMask) == input_result[1]


@pytest.mark.parametrize(
    "input_result",
    [
        (None, 0),
        ("all", 9),
        ("corner", 2),
        ("edge", 6),
        ("core", 1),
        (["corner", "edge"], 8),
        ("random4", 4),
        ("uniform+3", 3),
        (("core", "corner_random_2"), 3),
    ],
)
def test_storageMaskFunction_hexagonal(input_result):
    """Ensure valid storeStates are interpreted correctly."""
    config = THIS_DIR + "/data/hexagonal.yaml"
    S = Snowflake(N_vials=(3, 3, 1), storeStates=input_result[0], configPath=config)

    assert np.sum(S._storageMask) == input_result[1]


@pytest.mark.parametrize(
    "method_kwargs",
    [("to_frame", {}), ("sigmaCounter", {"time": 1})],
)
def test_mustRunFirst(method_kwargs):
    S = Snowflake()

    with pytest.raises(ValueError):
        getattr(S, method_kwargs[0])(**method_kwargs[1])


@pytest.fixture(scope="module")
def fakeS():
    """Fixture to share a fake Snowflake across tests.

    Fake in that states are manually provided to avoid rng issues.
    """
    config = THIS_DIR + "/data/square.yaml"
    S = Snowflake(N_vials=(2, 2, 1), storeStates="all", configPath=config)
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


@pytest.fixture(scope="module")
def fakeS_hex():
    """Fixture to share a fake Snowflake across tests.

    Fake in that states are manually provided to avoid rng issues.
    """
    config = THIS_DIR + "/data/hexagonal.yaml"
    S = Snowflake(N_vials=(2, 2, 1), storeStates="all", configPath=config)
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


def test_getvialgroup(S_331_all, S_331_all_hex, S_333_all, S_333_all_hex):
    """Test that vialgrouping masks are correctly computed."""
    assert np.where(S_331_all.getVialGroup("core"))[0] == 4
    assert all(
        np.where(S_331_all.getVialGroup(["core", "corner"]))[0] == [0, 2, 4, 6, 8]
    )

    assert S_331_all.getVialGroup("side").sum() == S_331_all.getVialGroup("edge").sum()
    assert (
        S_331_all_hex.getVialGroup("side").sum()
        == S_331_all_hex.getVialGroup("edge").sum()
    )

    with pytest.raises(ValueError):
        S_331_all.getVialGroup("rubbish")

    assert S_333_all.getVialGroup("core").sum() == 1
    assert S_333_all.getVialGroup("corner").sum() == 8
    assert S_333_all.getVialGroup("side").sum() == 6
    assert S_333_all.getVialGroup("edge").sum() == 12
    assert S_333_all.getVialGroup("all").sum() == 27

    assert S_333_all_hex.getVialGroup("core").sum() == 1
    assert S_333_all_hex.getVialGroup("corner").sum() == 4
    assert S_333_all_hex.getVialGroup("side").sum() == 22
    assert S_333_all_hex.getVialGroup("edge").sum() == 22
    assert S_333_all_hex.getVialGroup("all").sum() == 27


@pytest.mark.slow
def test_to_frame(fakeS, S_331_all):
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

    # making sure getVialGroup and to_frame have same
    # interpretation of groups
    df = S_331_all.to_frame()[0]
    groupVial = (
        df[["group", "vial"]].drop_duplicates().groupby("group")["vial"].unique()
    )

    for g, vials in groupVial.items():
        assert (np.where(S_331_all.getVialGroup(g))[0] == vials).all()


@pytest.mark.slow
def test_to_frame_hex(fakeS_hex, S_331_all_hex):
    """Test that conversion to dataframe works as intended."""
    stats_df, traj_df = fakeS_hex.to_frame(n_timeSteps=4)

    assert all(
        [col in stats_df.columns for col in ["group", "vial", "variable", "value"]]
    )
    assert stats_df.shape[0] == 12
    assert S_331_all_hex.getVialGroup("corner").sum() == 2

    assert (
        traj_df.loc[
            (traj_df.vial == 0) & (traj_df.state == "sigma") & (traj_df.Time == 3),
            "value",
        ].item()
        == 0.95
    )

    # making sure getVialGroup and to_frame have same
    # interpretation of groups
    df = S_331_all_hex.to_frame()[0]
    groupVial = (
        df[["group", "vial"]].drop_duplicates().groupby("group")["vial"].unique()
    )

    for g, vials in groupVial.items():
        assert (np.where(S_331_all_hex.getVialGroup(g))[0] == vials).all()


@pytest.fixture(scope="module")
def S_333_all():
    """Fixture to share Snowflake across tests."""
    config = THIS_DIR + "/data/square.yaml"
    S = Snowflake(N_vials=(3, 3, 3), configPath=config)
    return S


@pytest.fixture(scope="module")
def S_333_all_hex():
    """Fixture to share Snowflake across tests."""
    config = THIS_DIR + "/data/hexagonal.yaml"
    S = Snowflake(N_vials=(3, 3, 3), configPath=config)
    return S


@pytest.fixture(scope="module")
def S_331_all():
    """Fixture to share Snowflake across tests."""
    config = THIS_DIR + "/data/square.yaml"
    S = Snowflake(
        N_vials=(3, 3, 1),
        storeStates="all",
        initialStates=dict(temp=21),
        configPath=config,
    )
    S.run()
    return S


@pytest.fixture(scope="module")
def S_331_all_hex():
    """Fixture to share Snowflake across tests."""
    config = THIS_DIR + "/data/hexagonal.yaml"
    S = Snowflake(
        N_vials=(3, 3, 1),
        storeStates="all",
        initialStates=dict(temp=21),
        configPath=config,
    )
    S.run()
    return S


@pytest.fixture(scope="module")
def S_331_edge():
    config = THIS_DIR + "/data/square.yaml"
    """Fixture to share Snowflake across tests."""
    S = Snowflake(
        N_vials=(3, 3, 1),
        storeStates="edge",
        initialStates=dict(temp=21),
        configPath=config,
    )
    S.run()
    return S


@pytest.fixture(scope="module")
def S_331_edge_hex():
    config = THIS_DIR + "/data/hexagonal.yaml"
    """Fixture to share Snowflake across tests."""
    S = Snowflake(
        N_vials=(3, 3, 1),
        storeStates="edge",
        initialStates=dict(temp=21),
        configPath=config,
    )
    S.run()
    return S


def test_initialTemp(S_331_all, S_331_all_hex, S_331_edge, S_331_edge_hex):
    # we had defined initial temperature as 21
    assert all(S_331_all.X_T[:, 0] == 21)
    assert all(S_331_edge.X_T[:, 0] == 21)
    assert all(S_331_all_hex.X_T[:, 0] == 21)
    assert all(S_331_edge_hex.X_T[:, 0] == 21)


@pytest.fixture(scope="module")
def S_cnt():
    """Fixture with controlled nucleation that fails if turned on at -5."""
    d = {"int": 0, "ext": 0, "s0": 20, "s_sigma_rel": 0}
    S_cnt = Snowflake(
        storeStates="all",
        N_vials=(7, 7, 1),
        k=d,
        opcond=OperatingConditions(
            t_tot=4e4,
            cooling={"rate": 0.5 / 60, "start": 20, "end": -50},
            holding=[dict(duration=60, temp=-5), dict(duration=90 * 60, temp=-6)],
            cnTemp=None,
        ),
    )

    return S_cnt


@pytest.mark.slow
def test_cntWarning(S_cnt, capsys):
    """Test that warning appears if cnt fails (and doesn't otherwise)."""
    S_cnt.run()
    captured, _ = capsys.readouterr()

    assert "WARNING" not in captured

    # turn on cnt
    S_cnt.opcond.cnTemp = -5
    S_cnt.run()
    captured, _ = capsys.readouterr()

    # there a significant number of vials was not supercooled
    assert "WARNING" in captured


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
    config = THIS_DIR + "/data/square.yaml"
    S = Snowflake(N_vials=(3, 3, 1), configPath=config)
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

    # same checks but for the hexagonal setup
    config = THIS_DIR + "/data/hexagonal.yaml"
    S = Snowflake(N_vials=(3, 3, 1), configPath=config)
    mint, mext = S._buildInteractionMatrices()

    assert mint.shape == (9, 9)
    assert mint.nnz == 41
    assert mint[3, 0] == 1
    assert mint[0, 3] == 1
    assert mint[1, 1] == -4
    assert mint[3, 3] == -5
    assert mint[4, 4] == -6
    assert (np.sum(mint, axis=0) == 0).all()
    assert (np.sum(mint, axis=1) == 0).all()
