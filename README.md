<img src="https://github.com/SPL-ethz/snow/raw/main/docs_src/_media/snowLogo_v1.png" width="400">

![tag](https://img.shields.io/github/v/tag/SPL-ethz/snow) ![pypi](https://img.shields.io/pypi/v/ethz-snow) ![build](https://img.shields.io/github/workflow/status/SPL-ethz/snow/Python%20package) [![Coverage Status](https://coveralls.io/repos/github/SPL-ethz/snow/badge.svg?branch=feature/coveralls)](https://coveralls.io/github/SPL-ethz/snow?branch=feature/coveralls) ![](https://img.shields.io/github/stars/SPL-ethz/snow?style=social) ![](https://img.shields.io/github/watchers/SPL-ethz/snow?style=social) ![](https://img.shields.io/github/commit-activity/m/SPL-ethz/snow) ![](https://img.shields.io/github/issues-raw/SPL-ethz/snow) ![](https://img.shields.io/pypi/l/ethz-snow) ![](https://img.shields.io/badge/code%20style-black-black)

SNOW (**S**tochastic **N**ucleation **O**f **W**ater) is an open source modeling framework to simulate the freezing process in a large number of vials, tailored towards pharmaceutical freeze-drying processes. SNOW was developed as part of a research collaboration between ETH Zurich's Separation Processes Laboratory and Janssen (pharmaceutical companies of J&J). It is brought to you by Leif-Thore Deck and Dr. David Ochsenbein. 

## Description

SNOW is a model capable of simulating the entire freezing process for a batch of vials, from the liquid state to the completely frozen state. It combines balance equations, thermodynamics and a stochastic approach to ice nucleation to provide a comprehensive view of freezing and is entirely derived from first principles. 

In addition to the simulation of the thermal evolution of vials during freezing, SNOW tracks a number of characteristic quantities of the freezing process for each vial. These quantities are the nucleation time, the nucleation temperature and the solidification time. Both nucleation and solidification are correlated with attributes of the frozen product such as the mean crystal size and the product activity of an API. Variability in these quantities among vials may thus be used as surrogate for batch heterogeneity. SNOW consequently may be used in process design and optimization to identify suitable cooling protocols with sufficiently small batch heterogeneity. 

A more detailed description of the modeling framework is presented in a recent publication, which is publicly accessible under https://doi.org/10.1016/j.ijpharm.2021.121276 While version 1.0 is tailored towards simulations of the freezing stage in freeze-drying, version 1.1. was developed for pallet freezing; pallet freezing is for example used in the manufacturing of the Janssen COVID-19 vaccine, which served as case study for the model. Version 1.1 and the case study on the COVID-vaccine are discussed in detail in a recent pre-print at ChemRxiv: https://doi.org/10.26434/chemrxiv-2022-gnwhf

### Additional features currently supported (as of version 1.1)
- Simulation of systems comprising vials arranged in three spatial dimensions
- Arbitrary initial temperature of the vials
- Improvements in the numerical implementation (Second method to compare the initial amount of ice formed, faster data saving)

### Features supported as of version 1.0
- Simulation of the thermal evolution and ice formation for a batch with arbitrary number of vials
- Arbitrary cooling protocols (i.e., user may choose cooling rate, integrate holding steps and controlled nucleation)
- Tracking of nucleation times, temperatures and solidification times for all vials
- Stochastic nucleation in the form of a Monte Carlo approach as well as controlled nucleation in the form of forced initiation of nucleation at a certain point in time
- Cubic geometry of vial and rectangular arrangement on the shelf

### Features in preparation
- Hexagonal arrangement of vials on the shelf for freeze-drying
- Additional modes of heat transfer (e.g. thermal radiation)
- Alternative vial / vessel geometries
- Spatial simulation of freezing within a single vessel

## Bug reports
Please report bugs as [Github issues](https://github.com/SPL-ethz/snow/issues/new?assignees=ltdeck&labels=bug) or via [Email](mailto:deckl@ethz.ch), preferably
including a screenshot that illustrates the problem.

Copyright (c) 2021-2022 Leif-Thore Deck, David Ochsenbein

<sub>The snow package logo has been designed using resources from Flaticon.com.</sub>
