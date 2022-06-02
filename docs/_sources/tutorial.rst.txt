========
Tutorial
========

This tutorial provides some basic explanations to help users getting started with the ethz_snow package. 

Installation
============

After downloading the package, you need to install it. (The requirements for the python version and dependent packages are listed in the setup.cfg file and will be read automatically by pip.)

Either, download/clone any stable or unstable version from `GitHub <https://github.com/SPL-ethz/snow/>` Then navigate to the folder and install the package via ``pip install -e .``

Or install the latest release version via ``pip install ethz_snow``.

First steps: General information 
================================

The package contains several files, two of which enable the simulation of freezing processes. **snowflake.py** is the basis version of the model, capable of simulating the freezing process of a shelf comprising an arbitrary number of vials exactly once. **snowfall.py** is setup to run such simulation repetitively, namely for **Nrep** times. The idea behind the naming, of course, is that Snowfall in reality comprises a large number of Snowflakes, all of which are unique. In the same way, every freezing process will be unique as a consequence of the stochasticity of primary ice nucleation. 

Both models require a number of operatingConditions and constants as input parameters, which may be adjusted by the user. We define operating conditions as those input parameters, that are related to the numerics of the model (e.g. time step, number of repetitions) or to technical aspects (e.g. heat transfer parameters, cooling protocol). Constants, on the other hand, refer to formulation specific parameters, such as the volume of a vial, its composition and the required individual physicochemical parameters of the individual components (including ice nucleation kinetics). 

Operating conditions are to be provided directly when calling the **snowflake** / **Snowfall** functions, while we specify the constants separately in a .yaml file. If no values are specified, the default configurations are used; these are stored in the **operatingConditions.py** and the **snowConfig_default.yaml** (see the next subsection for more details on the latter).

Constants
---------
Below is a copy of the yaml file used to configure the constants used in **Snowflake** / **Snowfall**.
If you wish to modify any (subset) of these constants, simply create a new yaml that follows the same schema in an arbitrary location and link to it via the **configPath** input of Snowflake.
A partial config is allowed (e.g., one only containing a new water:cp_w entry). The remaining keys will be taken from the default.

.. code-block:: yaml

   vial:
    geometry:
    # base shape of vial (currently only cube is accepted)
    shape: cube
    # length [m]
    length: 0.01
    # width [m]
    width: 0.01
    # height [m]
    height: 0.01
 water:
  # J/[Kkg] heat capacity liquid water
  cp_w: 4187
  # J/[Kkg] heat capacity ice 0°C
  cp_i: 2108
  # J/kg heat of fusion water
  Dh: 333550
 solution:
  # kg/m^3 density of water / assumed constant for all phases
  rho_l: 1000
  # K kg / mol cryoscopic constant
  k_f: 1.853
  # J/Kkg, heat capacity sucrose
  cp_s: 1240
  # molar mass sucrose [kg/mol]
  M_s: 0.3423
  # mass fraction of solute in solution
  solid_fraction: 0.05
  # °C equilibrium freezing temperature pure water
  T_eq: 0
  # initial temperature of the solution
  T_init: 20
 kinetics:
  kb: 1e-9 # m−3 s−1 K−b
  b: 12

Example
========

Let us consider we would like to run a simulation following the default parameters in the snowConfig_default.yaml, however with a specific set of operating conditions. Indeed, we want to study the freezing of a system with a slow, but variable shelf-to-vial heat transfer, 7 times 7 vials on the shelf, a cooling rate of 1 K/min and two holding steps at -5°C and -10°C. This relates to typical conditions for freezing in a lab freeze-dryer. 

We first import both Snowfall and the OperatingConditions:

.. code-block:: python

    from ethz_snow.snowfall import Snowfall
    from ethz_snow.operatingConditions import OperatingConditions

Then, we define the values of the four heat transfer coefficients: "int" refers to thermal interaction among vials, "ext" to thermal interaction of the edge vials with the environment. We neglect both effects and set the values to zero. "s0" refers to the mean shelf-to-vial heat transfer coefficient, "s_sigma_rel" refers to its relative variability. We assume a rather low value of 20 for "s0" and 0.1 for its variability. Note that the pre-defined unit of the heat transfer coefficients is W/m^2K and that there is currently no feature to change that unit.

.. code-block:: python

    d = {"int": 20, "ext": 0, "s0": 20, "s_sigma_rel": 0.1}

Next, we define the cooling protocol. Note that holding steps are defined separately. In terms of cooling, we set the cooling rate, the start temperature and the end temperature as follows. Here, the temperatures are defined in °C and the cooling rate in K/s. Typical values are in the range between 0.1 - 1.0 K/min. The start temperature typically is set to ambient temperature, while the final temperature may depend on the technical capabilities of the freezing device.  

.. code-block:: python

    c = {"rate": 0.5 / 60, "start": 20, "end": -50}

For the holding steps, we need to define duration and temperature of each step separately. Let us say, that both steps at -5°C and at -10°C take 90 min. Again, it is important that the time is based in s, thus we need to multiply with 60. Note that the program will automatically adjust the sequence of the holding steps in the way that they are in the order of decreasing holding temperatures.

.. code-block:: python

    h = [dict(duration=90*60, temp=-10), dict(duration=90*60, temp=-5)]

Next, let us think about the total time of the simulation that is required; this depends on the cooling and holding parameters as well as on heat transfer and on some of the formulation constants. It is thus not automatically calculated, but needs to be set. We recommend to provide at least one hour more than is required for the shelf to reach the final temperature. One may use the Snowflake simulation to test if the set time is sufficient. Here, let us set t_tot = 3e4:

.. code-block:: python

    op = OperatingConditions(t_tot=3e4, cooling=c, holding=h)

In case, we are interested in controlled nucleation, we can add the argument cnTemp to trigger nucleation at the end of a holding step. By defining

.. code-block:: python

    op = OperatingConditions(t_tot=3e4, cooling=c, holding = h, cnTemp = -5 )

we trigger nucleation at the end of the holding step at -5°C. Note that in the current version, controlled nucleation may only be defined at the end of a holding step.

Finally, we may define the Snowfall class. We set the pool_size parameter to the number of available workers and Nrep to a statistically relevant number. To fully capture the effects of the stochasticity of ice nucleation in a quantitative manner, we recommend Nrep > 1000. For a qualitative view, we set Nrep = 50:

.. code-block:: python

    S = Snowfall(pool_size=8, k=d, Nrep=50, N_vials=(7,7,1), opcond=op)

We then start the simulation via **S.run()** and may check whether it completed via **S.simulationStatus()**. In case we are only interested in a single repetition, the **Snowflake** class may be used instead. Compared to **Snowfall**, **Snowflake** does not require Nrep or pool_size as input. However, it is able to store information on the thermal evolution of all vials, which is a feature that was removed for **Snowfall** to increase computational performance. 

Simulation output
=================

After running the simulation, several information are stored that characterize the freezing process. Importantly, these are the **solidificationTimes()**, **nucleationTimes()**, and **nucleationTemperatures()**. These are also grouped based on position, allowing to understand potential differences among center, edge and corner vials. 

We may use **S.plot(what="T_nucleation")** to immediately get an understanding of the nucleation temperatures, and similarly for the other quantities. The plot function is also capable of showing trajectories, in case **Snowflake** is used instead of **Snowfall**. In this case, 

.. code-block:: python

    S.plot(what="T_nucleation")

will show the evolution of the temperatures as well as the shelf, which is a very useful first information for understanding the freezing process as well as a sanity check of the simulation outcome. Note that the plotting of trajectories is slow at the moment because of the way seaborn calculates the shaded area (representing all the trajectories).

Version 1.1. Pallet freezing
============================

The main additional feature of version 1.1 is the capability to simulate the freezing of systems with vials arranged in three spatial dimensions, e.g. in pallets. These pallets may comprise tens of thousands of vials and are commonly frozen in cold storage rooms over the course of days. Pallet simulations are initiated in the same way as two dimensional arrangements; however, the number of vials in the vertical (z) direction is to be set to a value larger than one. For example, a system with 40x36x18 vials may be setup via


.. code-block:: python

    S = Snowfall(pool_size=8,k=d,Nrep=128,N_vials=(40,36,18),opcond=op,dt=5, initialStates = initial)
    
    
Note that due to the geometry applied, the heat transfer settings for a pallet configuration may be different than for freezing on a shelf. Specifically, no shelf-to-vial heat transfer may be present and the external, i.e. the storage temperature is most often constant. However, the storage temperature is colder than the initial temperature of the vials; this difference between initial temperature and storage temperature is considered via the new option **initialStates**. A sample configuration may be

.. code-block:: python

    d = {"int": 10, "ext": 10, "s0": 0, "s_sigma_rel": 0} 
    c = {"rate": 1e-16, "start": -8, "end": -8} # rate equals zero yields divide by zero error
    initial = {"temp": 20}
    op = OperatingConditions(t_tot=6e6,cooling= c, holding =dict(duration=6e6,temp=-8) )
    
In order to simulate a constant storage temperature, an arbitrarily small cooling rate may be defined in addition with a holding step. In this way, the temperature is set for the entire process duration to the storage temperature, which is -8°C in this case. Note that due to the large system size, typically longer process durations have to be simulated for pallets compared to shelf freezing. 