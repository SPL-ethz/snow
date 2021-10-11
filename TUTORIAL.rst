========
Tutorial
========

This tutorial provides some basic explanations to help users getting started with the ethz_snow package. 

========
Installation
========

After downloading the package, you need to install it. The requirements for the python version and dependent packages are listed in the setup.cfg file. You may install the package via ``pip install -e .``

========
First steps: General information 
========

The package contains several files, two of which enable the simulation of freezing processes. **snowflake.py** is the basis version of the model, capable of simulating the freezing process of a shelf comprising an arbitrary number of vials exactly once. **snowfall.py** is setup to run such simulation repetitively, namely for **Nrep** times. The idea behind the naming, of course, is that snowfall in reality comprises a large number of snowflakes, all of which are unique. In the same way, every freezing process will be unique as a consequence of the stochasticity of primary ice nucleation. 

Both models require a number of operatingConditions and constants as input parameters, which may be adjusted by the user. We define operating conditions as those input parameters, that are related to the numerics of the model (e.g. time step, number of repetitions) or to technical aspects (e.g. heat transfer parameters, cooling protocol). Constants, on the other hand, refer to formulation specific parameters, such as the volume of a vial, its composition and the required individual physicochemical parameters of the individual components (including ice nucleation kinetics). 

Operating conditions are to be provided directly when calling the **snowflake** / **snowfall** functions, while we specify the constants separately in a .yaml file. If no values are specified, the default configurations are used; these are stored in the **operatingConditions.py** and the **snowConfig_default.yaml**.

========
Example
========

Let us consider we would like to run a simulation following the default parameters in the snowConfig_default.yaml, however with a specific set of operating conditions. Indeed, we want to study the freezing of a system with a slow, but variable shelf-to-vial heat transfer, 7 times 7 vials on the shelf, a cooling rate of 1 K/min and two holding steps at -5°C and -10°C. This relates to typical conditions for freezing in a lab freeze-dryer. 

We first import both Snowfall and the OperatingConditions:

``from ethz_snow.snowfall import Snowfall``
``from ethz_snow.operatingConditions import OperatingConditions``

Then, we define the values of the four heat transfer coefficients: "int" refers to thermal interaction among vials, "ext" to thermal interaction of the edge vials with the environment. We neglect both effects and set the values to zero. "s0" refers to the mean shelf-to-vial heat transfer coefficient, "s_sigma_rel" refers to its relative variability. We assume a rather low value of 20 for "s0" and 0.1 for its variability. Note that the pre-defined unit of the heat transfer coefficients is W/m^2K and that there is currently no feature to change that unit.

``d = {"int": 20, "ext": 0, "s0": 20, "s_sigma_rel": 0.1}``

Next, we define the cooling protocol. Note that holding steps are defined separately. In terms of cooling, we set the cooling rate, the start temperature and the end temperature as follows. Here, the temperatures are defined in °C and the cooling rate in K/s. Typical values are in the range between 0.1 - 1.0 K/min. The start temperature typically is set to ambient temperature, while the final temperature may depend on the technical capabilities of the freezing device.  

``c = {"rate": 0.5 / 60, "start": 20, "end": -50}``

For the holding steps, we need to define duration and temperature of each step separately. Let us say, that both steps at -5°C and at -10°C take 90 min. Again, it is important that the time is based in s, thus we need to multiply with 60. Also, the program automatically adjusts the sequence of the holding steps in the way that they are in the order of decreasing holding temperatures.

``h = [dict(duration=90*60, temp=-5), dict(duration=90*60, temp=-10)]``

Next, let us think about the total time of the simulation that is required; this depends on the cooling and holding parameters as well as on heat transfer and on some of the formulation constants. It is thus not automatically calculated, but needs to be set. We recommend to provide at least one hour more than is required for the shelf to reach the final temperature. One may use the snowflake simulation to test if the set time is sufficient. Here, let us set t_tot = 3e4:

``op = OperatingConditions(t_tot=3e4,cooling= c,holding = h )``

In case, we are interested in controlled nucleation, we can add the argument cnTemp to trigger nucleation at the end of a holding step. By defining

``op = OperatingConditions(t_tot=3e4,cooling= c,holding = h, cnTemp = -5 )``

we trigger nucleation at the end of the holding step at -5°C. Note that in the current version, controlled nucleation may only be defined at the end of a holding step.

Finally, we may define the snowfall class. We set the pool_size parameter to the number of cores of the processor and the Nrep to a statistically relevant number. To fully capture the effects of the stochasticity of ice nucleation in a quantitative manner, we recommend Nrep > 1000. For a qualitative view, we set Nrep = 50:

``S = Snowfall(pool_size=8,k=d,Nrep=50,N_vials=(7,7,1),opcond=op)``

We then start the simulation via ``S.run()`` and may check whether it completed via ``S.simulationStatus()``. In case we are only interested in a single repetition, the **snowflake** class may be used instead. Compared to **snowfall**, **snowflake** does not require Nrep or pool_size as input. However, it is able to store information on the thermal evolution of all vials, which is a feature that was removed for **snowfall** to increase computational performance. 

========
Simulation output
========

After running the simulation, several information are stored that characterize the freezing process. Importantly, these are the **solidificationTimes()**, **nucleationTimes()**, and **nucleationTemperatures()**. These are also grouped based on position, allowing to understand potential differences among center, edge and corner vials. 

We may use ``S.plot(what="T_nucleation")`` to immediately get an understanding of the nucleation temperatures, and similarly for the other quantities. The plot function is also capable of showing trajectories, in case **snowflake** is used instead of **snowfall**. In this case, 

``S.plot(what="T_nucleation")``

will show the evolution of the temperatures as well as the shelf, which is a very useful first information for understanding the freezing process as well as a sanity check of the simulation outcome. 


