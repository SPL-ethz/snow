"""Implement unit tests for constants module."""

import pytest
from ethz_snow import constants as scon
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_nestedDictUpdate():
    """Test nested dictionary is updated correctly."""
    d = {"a": 1, "b": 2, "c": {"d": 3}}
    u = {"c": {"d": 4}}

    dnew = scon._nestedDictUpdate(d, u)

    assert dnew["c"]["d"] == 4


def test_warningWrongYamlSchema(capsys):
    """Test schema error is warned about."""
    config = scon._loadConfig(THIS_DIR + "/data/partiallyValidConfig.yaml")

    # this captures the stdout (print messages)
    captured, _ = capsys.readouterr()

    # there are keys that don't belong there
    assert "WARNING" in captured

    # the one valid key was updated
    assert config["vial"]["geometry"]["length"] == 1

    const = scon.calculateDerived(THIS_DIR + "/data/partiallyValidConfig.yaml")

    # ensure that right type
    assert isinstance(const, dict)

    # check that volume calc is correct
    assert const["V"] == (
        config["vial"]["geometry"]["length"]
        * config["vial"]["geometry"]["width"]
        * config["vial"]["geometry"]["height"]
    )


def test_notImplementedShape():
    """Test not implemented shapes raise error."""
    # make sure that non-cubic shape raises exception
    with pytest.raises(NotImplementedError):
        config = scon.calculateDerived(THIS_DIR + "/data/notImplementedShape.yaml")
