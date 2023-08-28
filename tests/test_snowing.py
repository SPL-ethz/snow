"""Implement unit tests for Snowing class."""
import os
import pytest
import numpy as np
import ethz_snow.utils as Utils
from ethz_snow.snowing import Snowing
from ethz_snow import constants as scon
from ethz_snow.operatingConditions import OperatingConditions

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_notImplementedDimensionality():
    """Test to make sure that not correctly specified model dimensionality
    raises an error."""
    with pytest.raises(NotImplementedError):
        config = scon.calculateDerived(
            THIS_DIR + "/data/notImplementedDimensionality.yaml"
        )


def test_notImplementedConfiguration():
    """Test to make sure that not correctly specified freezing configuration
    raises an error."""
    with pytest.raises(NotImplementedError):
        config = scon.calculateDerived(
            THIS_DIR + "/data/notImplementedConfiguration.yaml"
        )


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
    """Check that value error is raised if results requested before
    simulation is carried out."""
    S = Snowing()

    with pytest.raises(AssertionError):
        S.results

    with pytest.raises(AssertionError):
        S.iceMassFraction

    with pytest.raises(AssertionError):
        S.temp

    with pytest.raises(AssertionError):
        S.time

    with pytest.raises(AssertionError):
        S.shelfTemp


def test_totalTimeTooShort_for_nucleation():
    """Check that total time of simulation is long enough for
    nucleation to occur."""
    d = {"int": 0, "ext": 0, "s0": 20, "s_sigma_rel": 0}
    c = {"rate": 1 / 60, "start": 20, "end": 19}
    op = OperatingConditions(t_tot=200, cooling=c)

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/homogeneous_shelf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/spatial_1D_shelf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/spatial_2D_shelf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()


def test_totalTimeTooShort_for_complete_solidification():
    """Check that total time of simulation is long enough for
    complete solidification."""
    d = {"int": 0, "ext": 0, "s0": 20, "s_sigma_rel": 0}
    c = {"rate": 10 / 60, "start": 20, "end": -10}
    op = OperatingConditions(t_tot=3600, cooling=c, cnTemp=-2)

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/homogeneous_shelf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/spatial_1D_visf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()

    with pytest.raises(ValueError):
        config = THIS_DIR + "/data/spatial_2D_visf.yaml"
        S = Snowing(k=d, opcond=op, configPath=config)
        S.run()


def test_final_values_homogeneous():
    """Check that the final values of temperature and ice mass
    fraction make sense."""

    d = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}
    c = {"rate": 10 / 60, "start": 20, "end": -50}
    op = OperatingConditions(t_tot=10 * 3600, cooling=c)

    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(k=d, opcond=op, configPath=config)
    S.run()
    assert S.shelfTemp[-1] == -50
    assert np.round(S.temp[-1], 1) == -50
    assert np.round(S.time[-1], 1) == 10


def test_final_values_1D():
    """Check that the final values of temperature and ice mass
    fraction make sense."""

    d = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}
    c = {"rate": 10 / 60, "start": 20, "end": -50}
    op = OperatingConditions(t_tot=100 * 3600, cooling=c)

    config = THIS_DIR + "/data/spatial_1D_shelf.yaml"
    S = Snowing(k=d, opcond=op, configPath=config)
    S.run()
    assert S.shelfTemp[-1] == -50
    assert np.round(S.temp[-1], 1).min() == -50
    assert np.round(S.time[-1], 1) == 100


def test_final_values_2D():
    """Check that the final values of temperature and ice mass
    fraction make sense."""

    d = {"int": 0, "ext": 0, "s0": 50, "s_sigma_rel": 0}
    c = {"rate": 10 / 60, "start": 20, "end": -50}
    op = OperatingConditions(t_tot=1000 * 3600, cooling=c)

    config = THIS_DIR + "/data/spatial_2D_jacket.yaml"
    S = Snowing(k=d, opcond=op, configPath=config)
    S.run()
    assert S.shelfTemp[-1] == -50
    assert np.round(S.temp[-1], 1).min() == -50
    assert np.round(S.time[-1], 0) == 1000


def test_single_simulation_plot_cdf():
    """Check that plot_cdf doesn't work for single simulations."""
    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(configPath=config)
    with pytest.raises(NotImplementedError):
        S.run()
        S.plot_cdf()


def test_single_simulation_categorical_plot():
    """Check that plot doesn't work for single simulations."""
    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(configPath=config)
    with pytest.raises(NotImplementedError):
        S.run()
        S.plot()


def test_multiple_simulation_plot_evolution():
    """Check that plot_evolution doesn't work for multiple simulations."""
    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(configPath=config, Nrep=2)
    with pytest.raises(NotImplementedError):
        S.run()
        S.plot_evolution()


def test_incorrect_key_plot_evolution():
    """Check that plot_evolution doesn't work for incorrect key."""
    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(configPath=config)
    with pytest.raises(ValueError):
        S.run()
        S.plot_evolution(what="random")


def test_incorrect_key_plot_cdf():
    """Check that plot_cdf doesn't work for incorrect key."""
    config = THIS_DIR + "/data/homogeneous_shelf.yaml"
    S = Snowing(configPath=config, Nrep=10)
    with pytest.raises(ValueError):
        S.run()
        S.plot_cdf(what="random")


""" Unit tests for Utils class containing helper functins for Snowing. """


def test_vapour_pressure_liquid():
    """Test to make sure that vapour pressure above liquid water is computed
    correctly. For this the calculated values are compared to the experimental
    values reported in literature.
    """

    # compute vapour pressures at three different temperatures
    pressure = Utils.vapour_pressure_liquid(np.array([0, 10, 20, 30]) + 273.15)

    # check the values
    assert (np.round(pressure, -1) == np.array([610, 1230, 2340, 4250])).all()


def test_vapour_pressure_solid():
    """Test to make sure that vapour pressure above ice is computed correctly.
    For this the calculated values are compared to the experimental
    values reported in literature."""

    # compute vapour pressures at three different temperatures
    pressure = Utils.vapour_pressure_solid(np.array([0, -10, -20, -30]) + 273.15)

    # check the values
    assert (np.round(pressure, -1) == np.array([610, 260, 100, 40])).all()


def test_vapour_flux():
    """Test to make sure that vapour flux during vacuum is computed correctly.
    For this the calculated values are compared to the experimental
    values reported in literature."""

    # two different temperatures
    T_l_1 = 0 + 273.15
    T_l_2 = -10 + 273.15

    # check the values
    assert (
        np.round(
            np.array(
                [
                    Utils.vapour_flux(
                        0.015,
                        2.99e-26,
                        1.38e-23,
                        100,
                        Utils.vapour_pressure_solid(T_l_1),
                        T_l_1,
                        T_l_1,
                    ),
                    Utils.vapour_flux(
                        0.015,
                        2.99e-26,
                        1.38e-23,
                        100,
                        Utils.vapour_pressure_solid(T_l_2),
                        T_l_2,
                        T_l_2,
                    ),
                    Utils.vapour_flux(
                        0.004,
                        2.99e-26,
                        1.38e-23,
                        100,
                        Utils.vapour_pressure_solid(T_l_2),
                        T_l_2,
                        T_l_2,
                    ),
                ]
            ),
            4,
        )
        == np.array([87e-4, 28e-4, 7e-4])
    ).all()


def test_colormap_colors():
    """Make sure that the colormap gives out the correct starting color for a
    given spatial array z."""
    colors, my_cmap, norm = Utils.colormap(z=np.linspace(0, 1, 3))
    assert (colors[0, :] == np.array([0.0, 0.0, 0.0, 1.0])).all()
