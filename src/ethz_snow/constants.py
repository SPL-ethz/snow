"""Calculate derived constants and shortcut definitions."""

# Loading yaml config file
import pkg_resources
import yaml

from collections.abc import Mapping

from typing import Optional

# TODO:
# - config overwrite warning if key does not exist in default
# - type hints
# - doc strings
# - unit testing
# - ensure that this works in new installation


def _getAllKeys(dl, keys_list):
    if isinstance(dl, dict):
        keys_list += dl.keys()
        map(lambda x: _getAllKeys(x, keys_list), dl.values())
    elif isinstance(dl, list):
        map(lambda x: _getAllKeys(x, keys_list), dl)


def _nestedDictUpdate(d, u):
    # courtesy of stackoverflow (Alex Martelli)
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = _nestedDictUpdate(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def _loadConfig(fpath: Optional[str] = None) -> dict:
    defaultConfig_fpath = pkg_resources.resource_filename(
        "ethz_snow", "config/snowConfig_default.yaml"
    )

    with open(defaultConfig_fpath) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    if fpath is not None:
        with open(fpath) as f:
            customConfig = yaml.load(f, Loader=yaml.FullLoader)
        config = _nestedDictUpdate(config, customConfig)

    return config


def calculateDerived(fpath: Optional[str] = None):
    # the below code uses somewhat clunky casting
    # it's needed because pyyaml parses certain numbers
    # in scienfitic notation as strings (YAML 1.1 vs 1.2 I suppose)
    # I'm too lazy to write something more sophisticated.

    config = _loadConfig(fpath)

    # derived properties of vial
    if not config["vial"]["geometry"]["shape"].startswith("cub"):
        raise NotImplementedError(
            (
                f'Cannot handle shape "{config["vial"]["geometry"]["shape"]}". '
                + "Only cubic shape is supported at this moment."
            )
        )
    A = float(config["vial"]["geometry"]["a"]) * float(config["vial"]["geometry"]["b"])
    V = A * float(config["vial"]["geometry"]["c"])

    # derived properties of solution
    cp_s = float(config["solution"]["cp_s"])
    solid_fraction = float(config["solution"]["solid_fraction"])
    cp_w = float(config["water"]["cp_w"])
    cp_i = float(config["water"]["cp_i"])

    Dcp = cp_i - cp_w
    cp_solution = (
        solid_fraction * cp_s + (1 - solid_fraction) * cp_w
    )  # heat capacity of solution

    # Shortcut definitions
    mass = float(config["solution"]["rho_l"]) * V
    hl = mass * cp_solution
    depression = (
        float(config["solution"]["k_f"])
        / float(config["solution"]["M_s"])
        * (solid_fraction / (1 - solid_fraction))
    )

    alpha = (
        -mass * float(config["water"]["Dh"]) * (1 - solid_fraction)
    )  # used for sigma time step
    beta_solution = depression * mass * cp_solution

    return beta_solution
