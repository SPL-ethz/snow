=========
Changelog
=========

Version 1.1
===========

- Support for 3D configurations, typically of use to simulate the freezing of vials for storage in pallets
- Enhancements:
    - Initial Temp of vials can now be set independently (used to be defined by operating condition)
    - Fixing a minor bug where freezing point depression was not considered
    - Allow for cooling rate of 0
    - Additional warnings if holding and cooling times exceed the total time
    - Other, minor bug fixes

Version 1.0
===========

- Implementation of freezing model for the freezing stage in freeze-drying of a batch of vials on a shelf in python
- Currently supported features: 
    - Simulation of the freezing process, i.e. of thermal evolution and of ice formation, for a batch with arbitrary number of vials
    - Arbitrary cooling protocols (i.e., user may choose cooling rate, integrate holding steps and controlled nucleation)
    - Tracking of nucleation times, nucleation temperatures and solidification times for all vials
    - Stochastic nucleation in the form of a Monte Carlo approach as well as controlled nucleation in the form of forced initiation of nucleation at a certain point 
    - Cubic geometry of vial and rectangular arrangement on the shelf


