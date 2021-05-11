'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''

## Geometry and materials
mylambda = 0.8*0.05 # W/[mK] thermal conductivity of glass / times 0.1 to take into account that the actual contact area is smaller
thickness = 2e-3 # [m] thickness of glass (2 mm, usually 2 layers of 1 mm)
leng = 0.01 # [m] side length of cubic vial volume
A = leng*leng
V = A*leng

## Sample properties
rho_l = 1000 # kg/m^3 density of water / assumed constant for all phases
k_f = 1.853 # K kg / mol cryoscopic constant
cp_w = 4187 # J/[Kkg] heat capacity liquid water
cp_i = 2108 # J/[Kkg] heat capacity ice 0°C
cp_s = 1240 # J/Kkg, heat capacity sucrose
solid_fraction = 0.05 # mass fraction of solute in solution
Dcp = cp_i - cp_w
cp_solution = solid_fraction*cp_s + (1-solid_fraction)*cp_w # heat capacity of solution
T_init = 20 # °C initial temperature
T_eq = 0 #°C equilibrium freezing temperature
T_nuc = -10 # °C nucleation temperature
Dh = 333550 # J/kg heat of fusion water 

#Nucleation parameters (Braatz)
kb = 1e-9 # m−3 s−1 K−b value reduced compared to Braatz paper (he had a value of 10)
b = 12

#Sublimation parameters (sublimation is not part of the model anymore)
C_emp = 0.1e-6 # Empirical parameter for Kurz and Fisher model
t_factor = 0.225 # tortuosity factor tau^2/epsilon
Mw = 0.01802 # kg/mol molar mass of water
M_s = 342.3/1000 # molar mass sucrose
Rg = 8.3144 # gas constant

T_sub0 = 225 # Initial temperature of sublimating vials
T_shelf = -10 + 273.15
pwc = 13 #Chamber pressure in Pa
drho = rho_l*(1-solid_fraction) # density change during sublimation
Dhsub = 2.83e6 # latent heat of sublimation [J/kg]

pcrit = 611.657 # [Pa] Critical pressure
Tcrit = 273.16 # [K] Critical temperature

## Shortcut definitions
mass = rho_l*V
hl = mass *cp_solution
hs = mass *cp_i
depression = k_f/M_s *(solid_fraction/(1-solid_fraction))

sigma_equivalent = Dh/cp_w # temperature rise required for complete solidification
sigma_solid = sigma_equivalent *cp_w /cp_i  # similar temperature rise for a frozen sample, i.e. with solid heat capacity
