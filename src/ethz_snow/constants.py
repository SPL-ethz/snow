"""Calculate derived constants and shortcut definitions."""

# Loading yaml config file
import pkg_resources
import yaml
from numpy import pi

from collections.abc import Mapping

from typing import Optional, List, Any

from .__init__ import __citation__, __bibtex__


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
    # kb = float(config["kinetics"]["kb"]) Changed to vial dependent nuc model (11.03.2024)
    a = float(config["kinetics"]["a"])
    b = float(config["kinetics"]["b"])
    c = float(config["kinetics"]["c"])

    rho_l = float(config["solution"]["rho_l"])
    height = float(config["vial"]["geometry"]["height"])
    diameter = float(config["vial"]["geometry"]["diameter"])

    # snowing parameters
    dimensionality = str(config["snowing_parameters"]["dimensionality"])
    configuration = str(config["snowing_parameters"]["configuration"])

    # snowfall parameters
    vial_arrangement = str(config["snowfall_parameters"]["vial_arrangement"])

    # check configuration
    if configuration not in {"shelf", "VISF", "jacket"}:
        raise NotImplementedError(
            (
                f'Configuration "{configuration}" '
                + 'not correctly specified, use "shelf", "jacket" or "VISF".'
            )
        )

    # check vial arrangement
    if vial_arrangement not in {"hexagonal", "square"}:
        raise NotImplementedError(
            (
                f'Vial arrangment "{vial_arrangement}" '
                + 'not correctly specified, use "hexagonal" or "square".'
            )
        )

    # check model dimensionality
    if dimensionality not in {"homogeneous", "spatial_1D", "spatial_2D"}:
        raise NotImplementedError(
            (
                f'Model dimensionality "{dimensionality}" '
                + "not correctly specified, use homogeneous, spatial_1D or spatial_2D."
            )
        )

    # derived properties of vial for cubic shape
    if config["vial"]["geometry"]["shape"].startswith("cub"):
        A = float(config["vial"]["geometry"]["length"]) * float(
            config["vial"]["geometry"]["width"]
        )

    # vial properties for cylindrical shape, commented out for now
    # elif config["vial"]["geometry"]["shape"].startswith("cyl"):
    #     geoKeys = config["vial"]["geometry"].keys()
    #     if ("radius" in geoKeys) and (config["vial"]["geometry"]["radius"] is not None):
    #         # radius is specified (not in defaults!)
    #         cyl_r = float(config["vial"]["geometry"]["radius"])
    #     elif ("diameter" in geoKeys) and (
    #         config["vial"]["geometry"]["diameter"] is not None
    #     ):
    #         # diameter is specified
    #         cyl_r = float(config["vial"]["geometry"]["diameter"]) / 2
    #     elif ("length" in geoKeys) and ("width" in geoKeys):
    #         # length and width specified
    #         assert float(config["vial"]["geometry"]["length"]) == float(
    #             config["vial"]["geometry"]["width"]
    #         ), "Length and width must match for cylinders."

    #         # interpret as diameter
    #         cyl_r = float(config["vial"]["geometry"]["length"]) / 2
    #     elif "length" in geoKeys:
    #         # only length specified
    #         cyl_r = float(config["vial"]["geometry"]["length"]) / 2
    #     elif "width" in geoKeys:
    #         # only width specified
    #         cyl_r = float(config["vial"]["geometry"]["width"]) / 2

    #     A = pi * cyl_r**2

    else:
        raise NotImplementedError(
            (
                f'Cannot handle shape "{config["vial"]["geometry"]["shape"]}". '
                + 'Only "cubic" shape is supported. In Snowing, the shape is automatically rewritten to "cylinder".'
            )
        )

    # # vial arrangement is only considered in the homogeneous model
    # if dimensionality in {"homogeneous", "spatial_1D", "spatial_2D"}:
    #     print(
    #         (
    #             f'WARNING: Vial arrangement "{vial_arrangement}" '
    #             + "is relevant only for Snowfall and Snowflake simulations."
    #         )
    #     )

    # volume
    V = A * height

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
    mass = rho_l * V

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

    T_eq_l = T_eq - depression

    # latent heat of fusion
    Dh = float(config["water"]["Dh"])

    # needed for solving nucleation
    k_f = float(config["solution"]["k_f"])
    M_s = float(config["solution"]["M_s"])

    # general constants
    sigma_B = float(config["general"]["sigma_B"])
    k_B = float(config["general"]["k_B"])

    # mass of solute and water in the solution
    mass_solute = mass * solid_fraction
    mass_water = mass * (1 - solid_fraction)

    constVars = [
        "dimensionality",
        "vial_arrangement",
        "T_eq",
        "T_eq_l",
        "a",
        "b",
        "c",
        "A",
        "V",
        "cp_s",
        "rho_l",
        "solid_fraction",
        "cp_w",
        "cp_i",
        "cp_solution",
        "mass",
        "hl",
        "depression",
        "alpha",
        "beta_solution",
        "Dh",
        "k_f",
        "M_s",
        "sigma_B",
        "k_B",
        "mass_solute",
        "mass_water",
        "height",
        "diameter",
        "configuration",
    ]

    # check that spatial model is chosen for VISF
    if configuration == "VISF":
        if dimensionality == "homogeneous":
            raise NotImplementedError(
                (f'For simulating "{configuration}" ' + "a spatial model is required.")
            )

        # set up additional parameters for VISF
        p_vac = float(config["VISF"]["p_vac"])
        kappa = float(config["VISF"]["kappa"])
        Dh_evaporation = float(config["VISF"]["Dh_evaporation"])
        m_water = float(config["VISF"]["m_water"])
        t_vac_start = float(config["VISF"]["t_vac_start"])
        t_vac_duration = float(config["VISF"]["t_vac_duration"])

        constVars.extend(
            [
                "p_vac",
                "kappa",
                "Dh_evaporation",
                "m_water",
                "t_vac_start",
                "t_vac_duration",
            ]
        )

    # check that spatial model is chosen for jacket
    elif configuration == "jacket":
        if dimensionality != "spatial_2D":
            raise NotImplementedError(
                (f'For simulating "{configuration}" ' + "a 2D model is required.")
            )

        # set up additional parameters for jacket-ramped freezing
        air_gap = float(config["jacket"]["air_gap"])
        lambda_air = float(config["jacket"]["lambda_air"])

        constVars.extend(["lambda_air", "air_gap"])

    # warn that cylindrical gemoetry instead of cubic geometry is used for spatial model
    if dimensionality != "homogeneous":
        # print(
        #     (
        #         f'WARNING: For simulating a "{dimensionality}" '
        #         + "model a cylindrical geometry is required. Shape is automatically rewritten to cylinder."
        #     )
        # )

        # this is not a lumped capacitance model
        # effective heat conductivity (for non-homog. vial temps)
        lambda_s = float(config["solution"]["lambda_s"])
        lambda_w = float(config["water"]["lambda_w"])
        lambda_i = float(config["water"]["lambda_i"])

        # effective heat conductivity of solution
        lambda_solution = solid_fraction * lambda_s + (1 - solid_fraction) * lambda_w

        # extend the list of constant variables
        constVars.extend(
            [
                "lambda_s",
                "lambda_i",
                "lambda_w",
                "lambda_solution",
            ]
        )

    # bundle things into a dict now
    # don't do it earlier for readability
    for myvar in constVars:
        const[myvar] = locals()[myvar]

    return const
