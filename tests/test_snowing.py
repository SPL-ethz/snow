"""Implement unit tests for Snowing class."""
import os
import pytest
from ethz_snow.snowing import Snowing
from ethz_snow import constants as scon
from ethz_snow.operatingConditions import OperatingConditions

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_notImplementedDimensionality():
    """Test to make sure that not correctly specified model dimensionality raises an error."""
    with pytest.raises(NotImplementedError):
        config = scon.calculateDerived(
            THIS_DIR + "/data/notImplementedDimensionality.yaml"
        )


def test_notImplementedConfiguration():
    """Test to make sure that not correctly specified freezing configuration raises an error."""
    with pytest.raises(NotImplementedError):
        config = scon.calculateDerived(
            THIS_DIR + "/data/notImplementedConfiguration.yaml"
        )


def test_IncorrectlySpecifiedInputDimensionality():
    """Test to make sure that incorrectly specified Snowing inputs raise error."""
    with pytest.raises(NotImplementedError):
        S = Snowing(dimensionality="random")


def test_IncorrectlySpecifiedInputConfiguration():
    """Test to make sure that incorrectly specified Snowing inputs raise error."""
    with pytest.raises(NotImplementedError):
        S = Snowing(configuration="random")


def test_k_must_be_fine():
    """Test that k is checked properly."""
    with pytest.raises(TypeError):
        S = Snowing(k=[1, 2, 3])

    # make sure k is checked w.r.t. its keys
    with pytest.raises(ValueError):
        S = Snowing(k={"int": 1, "ext": 1, "s0f": 1})

    # in Snowing only s0 is required, make sure it raises error if not there
    with pytest.raises(ValueError):
        S = Snowing(k=dict(int=1, ext=1))


def test_opcond_must_be_operatingcondition():
    """Check that opcond type is verified."""
    with pytest.raises(TypeError):
        S = Snowing(opcond=dict(a=1))


def test_mustRunFirst():
    """Check that value error is raised if results requested before simulation is carried out."""
    S = Snowing()

    with pytest.raises(AssertionError):
        S.getResults


def test_totalTimeTooShort():
    """Check that total time of simulation is long enough for nucleation to occur."""
    d = {"int": 0, "ext": 0, "s0": 20, "s_sigma_rel": 0}
    c = {"rate": 1 / 60, "start": 20, "end": 19}
    h = []
    op = OperatingConditions(t_tot=10, cooling=c)
    S = Snowing(k=d, opcond=op)

    with pytest.raises(ValueError):
        S.run()
