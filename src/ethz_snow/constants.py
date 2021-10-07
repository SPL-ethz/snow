"""Calculate derived constants and shortcut definitions."""

# Loading yaml config file
import pkg_resources
import yaml

from collections.abc import Mapping

from typing import Optional, List, Any

# TODO:
# - unit testing for this capability
# - ensure that this works in a from-scratch installation
# - documentation for users on which config options exist and how to do


def __getAllKeys_gen(dl: Any) -> list:
    """Help recursive key search with generator.

    Args:
        dl (Any): Some dictionary or dict value.

    Returns:
        list: The list of keys in dl.

    Yields:
        Iterator[list]: A generator of keys.
    """
    if isinstance(dl, dict):
        for val in dl.values():
            yield from __getAllKeys_gen(val)

        yield list(dl.keys())


def _getAllKeys(dl: dict) -> List[str]:
    """Return all keys in a potentially nested dict.

    Args:
        dl (dict): A dictionary of dictionaries of arbitrary depth.

    Returns:
        List[str]: All keys.
    """
    keys_list = list(__getAllKeys_gen(dl))

    return [key for subl in keys_list for key in subl]


def _nestedDictUpdate(d: dict, u: dict) -> dict:
    """Update a nested dictionary with another.

    Both dictionaries come from the config yaml.
    We only want to update the entries in the custom
    config (u) that are different from the default (d).
    Args:
        d (dict): The reference nested dictionary.
        u (dict): The nested dictionary containing the updates.

    Returns:
        dict: The updated nested dictionary.
    """
    # courtesy of stackoverflow (Alex Martelli)
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = _nestedDictUpdate(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def _loadConfig(fpath: Optional[str] = None) -> dict:
    """Load the default config file and the custom one.

    Args:
        fpath (Optional[str], optional): The filepath of the custom
            config file. Defaults to None.

    Returns:
        dict: The loaded config as dict.
    """
    # the default config is listed as part of the package data
    # therefore we can access it with pkg_resources
    defaultConfig_fpath = pkg_resources.resource_filename(
        "ethz_snow", "config/snowConfig_default.yaml"
    )

    with open(defaultConfig_fpath) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    if fpath is not None:

        with open(fpath) as f:
            customConfig = yaml.load(f, Loader=yaml.FullLoader)

        # ensure that custom updates are
        # subset of the valid set of keys
        # (there is an edge case where the user
        # uses a valid key nested in the wrong place
        # we currently do not control for that)
        validKeys = _getAllKeys(config)
        newKeys = _getAllKeys(customConfig)

        diffSet = set(newKeys) - set(validKeys)
        if len(diffSet) > 0:
            print(
                (
                    f"WARNING: Custom config key{'s'*(len(diffSet) > 1)} {diffSet} "
                    + f"{'is' if (len(diffSet) == 1) else 'are'} not valid "
                    + "and will be ignored."
                )
            )

        config = _nestedDictUpdate(config, customConfig)

    return config


def calculateDerived(fpath: Optional[str] = None) -> dict:
    """Compute the constants needed for Snowflake. Derive where needed.

    Args:
        fpath (Optional[str], optional): The filepath of the custom
            config file. Defaults to None.

    Raises:
        NotImplementedError: Vial geometry is not cubic.

    Returns:
        dict: A dictionary of constants.
    """
    # the below code uses somewhat clunky casting
    # it's needed because pyyaml parses certain numbers
    # in scienfitic notation as strings (YAML 1.1 vs 1.2 I suppose)
    # I'm too lazy to write something more sophisticated.

    config = _loadConfig(fpath)

    const = dict()
    # copy directly from yaml
    T_eq = float(config["solution"]["T_eq"])
    kb = float(config["kinetics"]["kb"])
    b = float(config["kinetics"]["b"])

    # derived properties of vial
    if not config["vial"]["geometry"]["shape"].startswith("cub"):
        raise NotImplementedError(
            (
                f'Cannot handle shape "{config["vial"]["geometry"]["shape"]}". '
                + "Only cubic shape is supported at this moment."
            )
        )
    A = float(config["vial"]["geometry"]["length"]) * float(
        config["vial"]["geometry"]["width"]
    )
    V = A * float(config["vial"]["geometry"]["height"])

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

    # bundle things into a dict now
    # don't do it earlier for readability
    constVars = [
        "T_eq",
        "kb",
        "b",
        "A",
        "V",
        "cp_s",
        "solid_fraction",
        "cp_w",
        "cp_i",
        "cp_solution",
        "mass",
        "hl",
        "depression",
        "alpha",
        "beta_solution",
    ]

    for myvar in constVars:
        const[myvar] = locals()[myvar]

    return const
