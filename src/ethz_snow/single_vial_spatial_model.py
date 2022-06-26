#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 10:49:01 2021

@author: andrazkosir
"""

"""
MODELLING THE FREEZING STAGE IN A FREEZE-DRYING PROCESS.

Author: Andraz Kosir
Date: 10.10.2021

Aim:
    
Simulating the freezing stage in a 1D problem for a
single vial on a cooling shelf (only the heat flux from the shelf  and the
heat flux from the surroundings to the top surface of the product inside the
vial is considered). The temperature of the mixture inside the vial is
a function of vertical position z in the vial. Glass temperature is not
modelled, instead the product is directly in contact with the shelf or the
surroundings.

Structure:
   
Code is separated in 3 parts: (i) cooling stage, (ii) ice
nucleation stage and (iii) solidification stage; The heat conduction
equations are solved with the Finite Difference Method (FDM) using the 
explicit Euler time discretization (FTCS scheme).

"""

# -------------------------------------------------------------------------- #
# PREAMBLE & REQUIRED MODULES
# -------------------------------------------------------------------------- #

""" Modules. """

import math
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import simps
from matplotlib.colors import LinearSegmentedColormap

""" Plot settings. """

# matplotlib LaTeX interpreter and font settings for plots
plt.rcParams.update({
    "text.usetex" : True,
    "font.size" : 15,
    "mathtext.fontset" : "custom",
    "mathtext.it" : "computer modern:italic",
    "mathtext.rm"  : "computer modern:normal",
    "mathtext.sf": "computer modern:normal",
    "mathtext.tt": "computer modern:normal",
    "legend.fontsize" : 12,
    "axes.titlesize" : 15 
})

# remaining plot attributes
plt.rcParams.update({
    "axes.xmargin" : 0.15, "axes.ymargin" : 0.15,
    "axes.axisbelow" : True,
    "axes.grid.which": "both",
    "xtick.direction" : "in", "ytick.direction" : "in",
    "xtick.bottom" : True, "xtick.top" : True,
    "ytick.left" : True, "ytick.right" : True,
    "lines.linewidth" : 3,
    "lines.markersize" : 6,
    "grid.color": "lightgrey",
    "grid.linestyle" : (5, (3, 3)),
    "grid.linewidth" : 1,
    "hatch.color": "black",
    "xtick.major.size" : 6, "xtick.major.width" : 1,
    "xtick.minor.size" : 3, "xtick.minor.width" : 0.5, 
    "ytick.major.size" : 6, "ytick.major.width" : 1,
    "ytick.minor.size" : 3, "ytick.minor.width" : 0.5,
    "axes.edgecolor" : "black"
})

# -------------------------------------------------------------------------- #
# INPUT PARAMETERS
# -------------------------------------------------------------------------- #

""" Geometry and compositions. """

# vial geometry
H = 0.01 # [m] height of mixture in the vial
d = 0.01 # [m] diameter of the vial
A = 0.25*math.pi*d**2 # [m2] surface area of the vial bottom
V = A*H # [m3] volume of mixture in the vial

# mixture composition and mass
rho = 1000 # [kg/m3] density for all phases and concentrations of solute
m_v = V*rho # [kg] mass of mixture in the vial
w_s = 0.05 # mass fraction of solute
m_s = w_s*m_v # [kg] mass of solute
m_w = (1 - w_s)*m_v # [kg] mass of water

""" Thermal properties. """

# heat capacities
cp_w = 4187 # [J/kgK] heat capacity of liquid water
cp_i = 2108 # [J/kgK] heat capacity of ice at 0 °C
cp_s = 1240 # [J/kgK] heat capacity of solute (e.g., sucrose)
cp_eff = cp_s*w_s + cp_w*(1 - w_s) # effective specific heat capacity
lambda_w = 333550 # [J/kg] latent heat of fusion

# heat conductivities
k_w = 0.598 # [W/mK] heat conductivity of water
k_i = 2.25 # [W/mK] heat conductivity of ice
k_s = 0.126 # [W/mK] heat conductivity of sucrose
k_eff = w_s*k_s + (1 - w_s)*k_w # [W/mK] effective heat conductivity

# heat diffusivity
alpha = k_eff/(cp_eff*rho) # [m2/s] cooling-stage thermal diffusivity

""" Heat transfer parameters. """
# overall heat transfer coefficient from the shelf (bottom of the vial)
K_shelf = 50 # [W/m2K]
h_g = 0 # [W/m2K] heat trasnfer coefficient at the top of the vial
T_g = 280 # [K] air temperature in the vial at the top
T_wall = 280 # [K] wall temperature

""" Radiation heat transfer parameters. """
# stefan-boltzmann cosntant
sigma_B = 5.67e-8 # [W/m2K4]
F_shelf = 0 # [-] view factor for shelf-top 
F_wall = 0 # [-] view factor for wall-top

# initial temperatures & info on BC
T0 = 20 + 273.15 # [K] initial temperature
T_shelf_0 = 20 + 273.15 # [K] initial temperature of the shelf
T_shelf_final = -50 + 273.15 # [K] final temperature of the shelf
cooling_rate = -1/60 # [K/min] cooling rate

""" Nucleation kinetics. """

# nucleation kinetics
k_b = 1e-9 # [1/m3] nucleation prefactor
b = 12 # [/] nucleation exponent

# function for definition of nucleation rate
def nucleation_rate(T_eq_0, T):
    # J > 0, only when supercooled
    if any(T_eq_0 - T_j > 0 for T_j in T):
        J = k_b*(T_eq_0 - T)**b*(T_eq_0 - T > 0)
    # otherwise, J = 0
    else:
        J = 0.0*T
    return J

""" Freezing point depression: Blagden's model. """

# freezing point depression
T_m = 0 + 273.15 # [K] melting point temperature of pure water
k_f = 1.853 # [K/mol] cryoscopic constant of water
M_s = 342.3/1000 # [kg/mol] molar mass of sucrose
b_s = (1/M_s)*w_s/(1 - w_s) # [mol/kg] initial molality of solution
T_eq_0 = T_m - k_f*b_s # [K] initial equilibrium temperature

""" Numerical solution check: compute total removed heat. """
# total heat added/removed (initialization)
Q_tot = 0

# -------------------------------------------------------------------------- #
# SPATIAL DISCRETIZATION & PREALLOCATION FOR SAVING RESULTS
# -------------------------------------------------------------------------- #

""" Spatial discretization. """

# spatial discretization for all stages
Nz = 30 # [/] number of spatial grid points
dz = H/Nz # [m] size of spatial step
z = np.linspace(0, H, Nz) # [m] z-position array

# duration of cooling and solidification stages
time_cool_exp = 2*3600 # [s] expected duration of cooling

""" Preallocation of arrays for saving results from the cooling stage."""

# temperatures & ice mass fractions for cooling stage
N_save_cool = 1000
cooling_temp_results = np.zeros((N_save_cool, Nz))
cooling_ice_results = np.zeros((N_save_cool, Nz))
cooling_shelf_results = np.zeros(N_save_cool)
cooling_time_results = np.zeros(N_save_cool)
i_save = 0

# -------------------------------------------------------------------------- #
# SOLVING MODEL EQUATIONS USING THE FINITE DIFFERNCE METHOD (FTCS scheme)
# -------------------------------------------------------------------------- #

""" 1. Cooling stage. """

# temporal discretization
dt = 0.4*dz**2/alpha # [s] size of time step (determined by the CFL cond.)
t_cool_exp = np.linspace(0, time_cool_exp, round(time_cool_exp/dt)) # [s]
Nt_exp = len(t_cool_exp) # [/] number of temporal grid points
sigma = alpha*dt/dz**2 # [/] just a FDM short-cut

# arrays of temperature
T_old = np.zeros(Nz) + T0 # [K] temperature profile at t = T_old
T_new = T_old # [K] temperature profile at t = t_(n+1)

# initial temperature profile
T_init = T_old

# shelf temperature and coooling rate
T_shelf_cool = T_shelf_0 + cooling_rate*t_cool_exp # [K] shelf temperatures
T_shelf_cool[T_shelf_cool < T_shelf_final] = T_shelf_final # end of cooling

# cooling stage progressbar
bfmt = "{desc}: {percentage:.2f} % |{bar}| {n:.0f}/{total_fmt} \
[time: {elapsed}, ETA: {remaining}]"
pbar_cool = tqdm(total=Nt_exp, desc="Cooling stage", colour="red",
                 dynamic_ncols=True, disable=False, bar_format=bfmt, 
                 position=0, leave=True)

# finite difference method for solving the heat conduction equation
for i in range(Nt_exp):
    
    # update progressbar
    pbar_cool.update(1)
   
    " Boundary condition at the bottom, z = 0. (BACKWARD DIFFERENCING) "
    # overall heat flux
    q_overall = K_shelf*(T_shelf_cool[i] - T_old[0])
    # boundary condition on the bottom of the vial: ghost grid point
    T_bottom = T_old[0] + q_overall*dz/k_eff
    
    " Boundary condition at the top, z = H. (FORWARD DIFFERENCING)  "
    # convection heat flux
    q_convection = h_g*(T_g - T_old[0])
    # radiation heat flux (from the shelf)
    q_radiation_shelf = sigma_B*F_shelf*(T_shelf_cool[i]**4 - T_old[Nz-1]**4)
    # radiation heat flux (from the wall)
    q_radiation_wall = sigma_B*F_wall*(T_wall**4 - T_old[Nz-1]**4)
    # boundary condition in the top of the vial
    T_top = T_old[Nz-1] + (q_convection + q_radiation_shelf +
                           q_radiation_wall)*dz/k_eff
    
    " Bottom point. "
    # updating the temperature at the bottom
    T_new[0] = T_old[0] + sigma*(T_old[1] - 2*T_old[0] + T_bottom)     
    
    " Centre points. "
    # updating the bulk temperatures
    T_new[1:Nz-1] = T_old[1:Nz-1] + sigma*(T_old[2:Nz] - 
                               2*T_old[1:Nz-1] + T_old[0:Nz-2])

    " Top point. "
    # updating the top temperature
    T_new[Nz-1] = T_old[Nz-1] + sigma*(T_top - 2*T_old[Nz-1] + T_old[Nz-2])
    
    " Processing the results of the time iteration. "
    # update temperatures
    T_old = T_new
    
    " Stochastic nucleation for the entire vial. "
    # nucleation rate
    J_z = nucleation_rate(T_eq_0, T_new)
    # frequency of nucleation events
    K_v = A*simps(J_z, z)
    # probability of nucleation
    P_nucl = K_v*dt
    
    " Position dependent stochastic nucleation. "
    # frequency of nucleation events
    K_v_array = A*dz*J_z
    # probability of nucleation
    P_nucl_array = K_v_array*dt
    
    " Check for local supercooling. "
    # check the range of LCS = LocalSupercooling
    LCS_i = np.where(T_new < T_eq_0, 1, 0)
    # the not so supercool part
    LCS_i_r = np.where(LCS_i == 1, 0, 1)

    " Saving results of the cooling stage. "
    # sparsely saving results
    if np.mod(i, np.ceil(Nt_exp/N_save_cool)) == 0:
        # saving temperature results
        cooling_temp_results[i_save, :] = T_new
        # saving shelf temperature results
        cooling_shelf_results[i_save] = T_shelf_cool[i]
        # saving sparse time array
        cooling_time_results[i_save] = dt*i
        # update saving index
        i_save = i_save + 1
                
    " Positional based Monte-Carlo approach. "
    # array of random numbers between 0 and 1
    P_rand_array = np.random.rand(Nz)
    if (P_nucl_array > P_rand_array).any():
        # position index of nucleation
        nucl_position_index = np.where(P_nucl_array > P_rand_array)[0]
        # position in the vial where nucleation occured (can be more than 1)
        nucl_position = z[nucl_position_index]
        # ending index of cooling stage
        Nt_cool_end, i_save_end = i, i_save - 1
        # time of nucleation
        nucleation_time = Nt_cool_end*dt/3600
        # kill the progressbar
        pbar_cool.update(Nt_exp-i-1); pbar_cool.close()
        break

# actual duration of cooling stage
time_cool_end = Nt_cool_end*dt
t_cool_end = np.linspace(0, time_cool_end + dt, round(Nt_cool_end*dt))

# final solution and plot for the end of the cooling stage
T_cooling = T_new

# ice mass fraction during the cooling stage
w_ice_cooling = np.zeros(Nt_cool_end + 1)

""" 2. Ice nucleation stage. """

# general parameters
T_nucl = T_new
t_nucl = Nt_cool_end*dt
t_nucl_h = t_nucl/3600

# initialization of equilibrium temperature and ice masses
T_eq = T_new
m_i_nucl = np.zeros(Nz)

# solving the quadratic equation (vectorized (A = 1)
B = - T_m - T_nucl - lambda_w*m_w/(cp_eff*m_v)
C = (lambda_w*m_w*T_m/(cp_eff*m_v) - 
      m_s*(k_f/M_s)*lambda_w/(cp_eff*m_v) + T_m*T_nucl)

# equilibrium temperatures determined from quadratic equation (A = 1)
T_eq_sol_1 = 0.5*(-B - np.sqrt(B**2 - 4*C))
T_eq_sol_2 = 0.5*(-B + np.sqrt(B**2 - 4*C))
T_eq = T_new*LCS_i_r + T_eq_sol_1*LCS_i

# computed mass of ice from the quadratic equation
m_i_nucl_sol_1 = m_w - m_s*(k_f/M_s) / (T_m - T_eq_sol_1)
m_i_nucl_sol_2 = m_w - m_s*(k_f/M_s) / (T_m - T_eq_sol_2)
m_i_nucl = np.zeros(Nz)*LCS_i_r + m_i_nucl_sol_1*LCS_i

# mass fraction after nucleation
w_i_nucl = m_i_nucl/m_v

# saving results
cooling_temp_results[i_save_end, :] = T_eq
cooling_ice_results[i_save_end, :] = w_i_nucl

# removing all zero-entries in the saved results (stochastic nucleation)
cooling_temp_results = cooling_temp_results[:i_save_end, :]
cooling_ice_results = cooling_ice_results[:i_save_end, :]

# temperatures & ice mass fractions for solidification stage
N_save_solid = 1000
solid_temp_results = np.zeros((N_save_solid, Nz))
solid_ice_results = np.zeros((N_save_solid, Nz))
solid_shelf_results = np.zeros(N_save_solid)
solid_time_results = np.zeros(N_save_solid)
i_save = 0

""" 3. Solidification stage. """

# [s] duration of further cooling for solidification
time_solid = 2*3600 - t_nucl

# limiting value of alpha for computing dt from CFL condition
alpha_max = k_i/(cp_i*rho)

# temporal discretization
dt = 0.4*dz**2/alpha_max # [s] size of time step (CFL determined)
t_solid = np.linspace(0, time_solid, round(time_solid/dt)) # [s] time array
Nt_solid = len(t_solid) # [/] number of temporal grid points
time_solid_results = np.linspace(0, time_solid, N_save_solid)

# arrays of temperature
T_old = T_eq # [K] temperature profile at t = T_old
T_new = T_old # [K] temperature profile at t = t_(n+1)
T_shelf_solid = T_shelf_cool[Nt_cool_end] + cooling_rate*t_solid # [K] shelf temp.
T_shelf_solid[T_shelf_solid < T_shelf_final] = T_shelf_final # end of cooling

# array for mass fraction of ice in solidification
w_i_new = w_i_nucl

# solidification stage progressbar
pbar_solid = tqdm(total=Nt_solid, desc="Solidification stage", colour="red",
                  dynamic_ncols=True, disable=False, bar_format=bfmt,
                  position=0, leave=True)

# FDM sceheme for solving the heat conduction equation
for i in range(Nt_solid):
        
    # update progressbar
    pbar_solid.update(1)
    
    " Thermal properties. "
    # heat capacity
    cp_eff = cp_s*w_s + cp_i*w_i_new + cp_w*(1 - w_s - w_i_new)
    # heat conductivity
    k_eff = k_i*w_i_new + k_w*(1 - w_i_new)
    # heat diffusivity
    alpha = k_eff/(cp_eff*rho)
    # nonlinear capacitance term
    beta = lambda_w*k_f*m_s/(M_s*rho*V*cp_eff)
    # total nonlinear capacitance term
    BETA = np.ones(Nz)*LCS_i_r + (1 + beta/(T_old - T_m)**2)*LCS_i
    
    " Boundary condition at the bottom, z = 0. (FORWARD DIFFERENCING) "
    # overall heat flux
    q_overall = K_shelf*(T_shelf_solid[i] - T_old[0])
    # boundary condition on the bottom of the vial: ghost grid point
    T_bottom = T_old[0] + q_overall*dz/k_eff[0]
    
    " Boundary condition at the top, z = H. (BACKWARD DIFFERENCING) "
    # convection heat flux
    q_convection = h_g*(T_g - T_old[0])
    # radiation heat flux (from the shelf)
    q_radiation_shelf = sigma_B*F_shelf*(T_shelf_solid[i]**4 - T_old[Nz-1]**4)
    # radiation heat flux (from the wall)
    q_radiation_wall = sigma_B*F_wall*(T_wall**4 - T_old[Nz-1]**4)
    # boundary condition in the top of the vial
    T_top = T_old[Nz-1] + (q_convection + q_radiation_shelf +
                           q_radiation_wall)*dz/k_eff[Nz-1]

    " Bottom point. "
    # update the point at the bottom
    T_new[0] = (T_old[0] + 
        (dt/(cp_eff[0]*rho))*((k_eff[1] - k_eff[0])*(T_old[1] - T_bottom)/(4*dz**2) +
        k_eff[0]*(T_old[1] - 2*T_old[0] + T_bottom)/dz**2)*BETA[0]**(-1))
    
    " Center points. "
    # update the bulk points
    T_new[1:Nz-1] = (T_old[1:Nz-1] +
        (dt/(cp_eff[1:Nz-1]*rho))*((k_eff[2:Nz] - k_eff[0:Nz-2])*(T_old[2:Nz] - T_old[0:Nz-2])/(4*dz**2) +
        k_eff[1:Nz-1]*(T_old[2:Nz] - 2*T_old[1:Nz-1] + T_old[0:Nz-2])/dz**2)*BETA[1:Nz-1]**(-1))

    " Top point. "
    # update the temperature at the top
    T_new[Nz-1] = (T_old[Nz-1] +
        (dt/(cp_eff[Nz-1]*rho))*((k_eff[Nz-1] - k_eff[Nz-2])*(T_top - T_old[Nz-2])/(4*dz**2) + 
        k_eff[Nz-1]*(T_top - 2*T_old[Nz-1] + T_old[Nz-2])/dz**2)*BETA[Nz-1]**(-1))

    " Processing the results I: temperatures. "
    # update temperatures
    T_old = T_new

    " Check for local supercooling. "
    # check the range of LCS
    LCS_i = np.where(T_new < T_eq_0, 1, 0)
    # the not so supercool region
    LCS_i_r = np.where(LCS_i, 0, 1)
    # position of the supercooling region limit
    position = next(iter(np.where(LCS_i == 0)[0]), -1)

    " Processing the results II: ice mass fractions. "
    # ice mass distribution
    m_ice = np.zeros(Nz)*LCS_i_r + (m_w - m_s*(k_f/M_s)/(T_m - T_new))*LCS_i
    # ice mass fraction distribution
    w_i_new = m_ice/m_v
    
    " Saving results of the cooling stage. "
    # sparse saving
    if np.mod(i, np.ceil(Nt_solid/N_save_solid)) == 0:
        # temperature results
        solid_temp_results[i_save, :] = T_new
        # ice mass fraction results
        solid_ice_results[i_save, :] = w_i_new
        # shelf temperature results
        solid_shelf_results[i_save] = T_shelf_solid[i]
        # saving sparse time array
        solid_time_results[i_save] = t_nucl + dt*i
        # updating saving index
        i_save = i_save + 1

    " Compute solidification time: 95 % of total ice formed. "
    w_i_total = simps(m_ice, z)/(m_v*H)
    # check for solidification time, then break loop
    if w_i_total >= 0.95*(1 - w_s):
        t_sol, Nt_sol = dt*i/60, i_save

# kill the progressbar
pbar_solid.close()

# -------------------------------------------------------------------------- #
# PROCESSING THE RESULTS
# -------------------------------------------------------------------------- #

""" Array concatenations for all three stages. """

# temperature distributions in the product
temp_results = np.concatenate((cooling_temp_results[:i_save_end-1,:], 
                        solid_temp_results[:i_save-1,:]), axis=0)

# ice mass fraction results
ice_results = np.concatenate((cooling_ice_results[:i_save_end-1,:],
                        solid_ice_results[:i_save-1,:]), axis=0)

# total time array
time_results = np.concatenate((cooling_time_results[:i_save_end-1],
                        solid_time_results[:i_save-1]), axis=0)/3600

# total shelf temperature array
shelf_results = np.concatenate((cooling_shelf_results[:i_save_end-1],
                        solid_shelf_results[:i_save-1]), axis=0)

# -------------------------------------------------------------------------- #
# PLOTTING THE SOLUTION
# -------------------------------------------------------------------------- #

# custom color map (black -> gray -> green)
my_cmap = LinearSegmentedColormap.from_list("my_BBG", ["black",
                                        "darkgray", "forestgreen"])
colors = my_cmap(np.linspace(0, 1, Nz))
N_lines = 11
cmap = my_cmap
colors = cmap(np.linspace(0, 1, Nz))
cmaplist = [cmap(i) for i in range(cmap.N)]
my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0, 1000*H, N_lines)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# temperature evolution
plt.figure(dpi = 500)
for i in range(0, Nz, round(Nz/11)):
    plt.plot(time_results, temp_results[:,i] - 273.15, color=colors[i])
plt.plot(time_results, shelf_results[:] - 273.15, "k--", linewidth = 1.5)
plt.xlabel(r"$\mathrm{time}, t\ [\mathrm{h}]$")
plt.ylabel(r"$\mathrm{temperature}, T\ [\mathrm{^{\circ}C}]$")
plt.plot([t_nucl_h, t_nucl_h], [-500, 500], "m--", linewidth = 1.5)
plt.text(0.15, -10, r"$T_{\mathrm{sh}}$", fontsize=20, color="black")
plt.text(0.42, -48, r"$t^{\mathrm{nuc}}$", fontsize=20, color="magenta")
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
cbar = plt.colorbar(sm)
cbar.set_label(r"$\mathrm{position}, z\ [\mathrm{mm}]$", rotation = 90)
plt.title(r"$\mathbf{Temperature\ evolution}$")
plt.grid(axis = "both", color = "gray", linestyle = "dashed", alpha = 0.5)
plt.xlim(-0.1*2, 1.1*2)
plt.ylim(-55, 25)
plt.show()

# ice mass fraction evolution
plt.figure(dpi = 500)
for i in range(0, Nz, round(Nz/11)):
    plt.plot(time_results, ice_results[:,i], color=colors[i])
plt.xlabel(r"$\mathrm{time}, t\ [\mathrm{h}]$")
plt.ylabel(r"$\mathrm{ice\ mass\ fraction}, w_{\mathrm{i}}\ [-]$")
plt.title(r"$\mathbf{Ice\ mass\ fraction\ evolution}$")
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
cbar = plt.colorbar(sm)
cbar.set_label(r"$\mathrm{position}, z\ [\mathrm{mm}]$", rotation = 90)
plt.grid(axis = "both", color = "gray", linestyle = "dashed", alpha = 0.5)
plt.plot([-5, 5], [0.95, 0.95], "b--", linewidth = 1.5)
plt.plot([t_nucl_h, t_nucl_h], [-5, 5], "m--", linewidth = 1.5)
plt.text(-0.05, 0.85, r"$w_{\mathrm{s}} = 0.05$", fontsize=20, color="blue")
plt.text(0.42, 0.1, r"$t^{\mathrm{nuc}}$", fontsize=20, color="magenta")
plt.xlim(-0.1*2, 1.1*2)
plt.ylim(-0.1, 1.1)
plt.show()


