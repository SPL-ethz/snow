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
