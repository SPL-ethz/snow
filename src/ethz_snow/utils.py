import numpy as np
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import colorsys

# -------------------------------------------------------------------------- #
# Vacuum-induced surface freezing helper functions
# -------------------------------------------------------------------------- #

# More information on the following consititutive equations used when
# simulating VISF can be found in the following literature:
#
# Marek, R.; Straub, J. Analysis of the evaporation coefficient and
# the condensation coefficient of water. International Journal of Heat
# and Mass Transfer 2001, 44, 39–53.
# https://doi.org/10.1016/S0017-9310(00)00086-7
#
# Murphy, D. M.; Koop, T. Review of the vapour pressures of ice and
# supercooled water for atmospheric applications. Quarterly Journal of
# the Royal Meteorological Society 2005, 131, 1539–1565.
# https://doi.org/10.1256/qj.04.94


# vapour pressure estimation (temperature in K, pressure in Pa)
def vapour_pressure_liquid(T_liq):
    """A function to compute the water vapour pressure above
    the liquid product in a vial based on the product surface
    temperature.

    Args:
        T_liq (np.ndarray): Liquid product temperature at the surface.

    Returns:
        np.ndarray: Water vapour pressure above the product surface.
    """
    # pressure in Pa
    p_liq = np.exp(
        54.842763
        - 6763.22 / T_liq
        - 4.210 * np.log(T_liq)
        + 0.000367 * T_liq
        + np.tanh(0.0415 * (T_liq - 218.8))
        * (53.878 - 1331.22 / T_liq - 9.44523 * np.log(T_liq) + 0.014025 * T_liq)
    )
    # return vapour pressure
    return p_liq


# vapour pressure estimation (temperature in K, pressure in Pa)
def vapour_pressure_solid(T_sol):
    """A function to compute the water vapour pressure above
    the frozen product in a vial based on the frozen product
    surface temperature.

    Args:
        T_sol (np.ndarray): Frozen product temperature at the surface.

    Returns:
        np.ndarray: Water vapour pressure above the product surface.
    """
    # pressure in Pa
    p_sol = np.exp(
        9.550426 - 5723.265 / T_sol + 3.53068 * np.log(T_sol) - 0.00728332 * T_sol
    )
    # return vapour pressure
    return p_sol


# vapour pressure estimation (temperature in K, pressure in Pa)
def vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v):
    """A function to compute the water vapour flux at the top surface,
    which is exposed to vacuum during VISF.

    Args:
        kappa (float): Evaporation efficiency.
        m_water (float): Mass of a water molecule.
        k_B (float): Boltzmann constant.
        p_vac (np.ndarray): Chamber vacuum pressure.
        p_vap (np.ndarray): Vapour pressure.
        T_l (np.ndarray): Product temperature.
        T_v (np.ndarray): Vapour temperature.

    Returns:
        np.ndarray: FLux of water vapour at the top surface of the frozen
        or liquid product
    """
    # calculate flux
    N_w = (
        (2 / (2 - kappa))
        * np.sqrt(m_water * kappa**2 / (2 * np.pi * k_B))
        * (p_vap / np.sqrt(T_l) - p_vac / np.sqrt(T_v))
    )
    # return vapour flux
    return N_w


# -------------------------------------------------------------------------- #
# Plotting helper functions
# -------------------------------------------------------------------------- #


def colormap(z):
    """A function to create a colormap used to plot the time evolutions of
    temperature and ice mass fractions with respect to vertical positions
    in the vial.

    Args:
        z (np.ndarray): Discretized spatial domain in the vertical direction.
    """

    def lighten_color(color, amount=0.5):
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    my_cmap = LinearSegmentedColormap.from_list(
        "my_BBG", ["black", lighten_color("lightslategray", 0.75)]
    )
    Nz = len(z)
    colors = my_cmap(np.linspace(0, 1, Nz))
    # number of lines to plot the evolutions
    N_lines = 11
    colors = my_cmap(np.linspace(0, 1, Nz))
    cmaplist = [my_cmap(i) for i in range(my_cmap.N)]
    my_cmap = LinearSegmentedColormap.from_list("Custom cmap", cmaplist, my_cmap.N)
    bounds = np.linspace(0, 1, N_lines)
    norm = BoundaryNorm(bounds, my_cmap.N)

    return colors, my_cmap, norm
