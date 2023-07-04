import numpy as np
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import colorsys

# -------------------------------------------------------------------------- #
# Vacuum-induced surface freezing helper functions
# -------------------------------------------------------------------------- #


# vapour pressure estimation (temperature in K, pressure in Pa)
def vapour_pressure_liquid(T_liq):
    """_summary_

    Args:
        T_liq (_type_): _description_

    Returns:
        _type_: _description_
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
    """_summary_

    Args:
        T_sol (_type_): _description_

    Returns:
        _type_: _description_
    """
    # pressure in Pa
    p_sol = np.exp(
        9.550426 - 5723.265 / T_sol + 3.53068 * np.log(T_sol) - 0.00728332 * T_sol
    )
    # return vapour pressure
    return p_sol


# vapour pressure estimation (temperature in K, pressure in Pa)
def vapour_flux(kappa, m_water, k_B, p_vac, p_vap, T_l, T_v):
    """_summary_

    Args:
        kappa (_type_): _description_
        m_water (_type_): _description_
        k_B (_type_): _description_
        p_vac (_type_): _description_
        p_vap (_type_): _description_
        T_l (_type_): _description_
        T_v (_type_): _description_

    Returns:
        _type_: _description_
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
    """_summary_

    Args:
        z (_type_): _description_
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
