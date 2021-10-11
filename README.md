# ethz_snow

SNOW (**S**tochastic **N**ucleation **O**f **W**ater) is an open source modeling framework to simulate the freezing process in a large number of vials, tailored towards pharmaceutical freeze-drying processes. SNOW was developed as part of a research collaboration between ETH Zurich's Separation Processes Laboratory and Janssen (pharmaceutical companies of J&J). It is brought to you by Leif-Thore Deck and Dr. David Ochsenbein.  

## Description

SNOW is a model capable of simulating the entire freezing process for a batch of vials, from the liquid state to the completely frozen state. It combines balance equations, thermodynamics and a stochastic approach to ice nucleation to provide a comprehensive view of freezing and is entirely derived from first principles. 

In addition to the simulation of the thermal evolution of vials during freezing, SNOW tracks a number of characteristic quantities of the freezing process for each vial. These quantities are the nucleation time, the nucleation temperature and the solidification time. Both nucleation and solidification are correlated with attributes of the frozen product such as the mean crystal size and the product activity of an API. Variability in these quantities among vials may thus be used as surrogate for batch heterogeneity. SNOW consequently may be used in process design and optimization to identify suitable cooling protocols with sufficiently small batch heterogeneity. 

A more detailed description of the modeling framework will be presented in a scientific publication that is currently in preparation. 

### Features currently supported (as of version 1.0)
- Simulation of the thermal evolution and ice formation for a batch with arbitrary number of vials
- Arbitrary cooling protocols (i.e., user may choose cooling rate, integrate holding steps and controlled nucleation)
- Tracking of nucleation times, temperatures and solidification times for all vials
- Stochastic nucleation in the form of a Monte Carlo approach as well as controlled nucleation in the form of forced initiation of nucleation at a certain point in time
- Cubic geometry of vial and rectangular arrangement on the shelf

### Features in preparation
- Hexagonal arrangement of vials on the shelf
- Flexible vial geometry

## Bug reports
Please report bugs as [Github issues](https://github.com/SPL-ethz/snow/issues/new?assignees=ltdeck&labels=bug) or as via [Email](mailto:deckl@ethz.ch), preferably
including a screenshot that illustrates the problem.

Copyright (c) 2021 Leif-Thore Deck, David Ochsenbein
