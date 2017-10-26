# SiMRiv (R package) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SiMRiv)](https://cran.r-project.org/package=SiMRiv) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) [![Downloads](http://cranlogs.r-pkg.org/badges/SiMRiv)](https://cran.r-project.org/package=SiMRiv) [![Travis-CI Build Status](https://travis-ci.org/miguel-porto/SiMRiv.svg?branch=devel)](https://travis-ci.org/miguel-porto/SiMRiv.svg?branch=devel)

**Individual-Based, Spatially-Explicit Simulation and Analysis of Multi-State Movements in River Networks, Homogeneous, and Heterogeneous Landscapes in R.**

Available on [CRAN](https://cran.r-project.org/web/packages/SiMRiv/index.html), can be directly installed from within R.

Provides functions to generate and analyze spatially-explicit simulations of multi-state movements in heterogeneous landscapes, based on "resistance" rasters.
Although originally conceived and designed to fill the gap of softwares simulating spatially-explicit trajectories of species constrained to linear habitats or dendritic ecological networks (e.g., river networks), the simulation algorithm is built to be highly flexible and can be applied to any (aquatic, semi-aquatic or terrestrial) organism, independently of the landscape in which it moves.
Thus, the user will be able to use the package to simulate movements either in homogeneous landscapes, heterogeneous landscapes (e.g. semi-aquatic animal in a riverscape),
or even in highly contrasted landscapes (e.g. fish in a river network). The algorithm and its input parameters are the same for all cases, so that results are comparable.

Simulated trajectories can then be used as null models to test e.g. for species site fidelity (Powell 2000) and other movement ecology hypotheses (Nathan et al. 2008)
or to build predictive, mechanistic movement models (Moorcroft and Lewis 2006), among other things. The package should thus be relevant to explore a broad spectrum of
ecological phenomena, such as those at the interface of animal behaviour, landscape, spatial and movement ecology, disease and invasive species spread, and population dynamics.

This is the first released experimental version; do test before using in production.

