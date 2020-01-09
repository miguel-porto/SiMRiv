# SiMRiv (R package) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SiMRiv)](https://cran.r-project.org/package=SiMRiv) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) [![Downloads](http://cranlogs.r-pkg.org/badges/SiMRiv)](https://cran.r-project.org/package=SiMRiv) [![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/SiMRiv)](https://cran.r-project.org/package=SiMRiv) [![Travis-CI Build Status](https://travis-ci.org/miguel-porto/SiMRiv.svg?branch=devel)](https://travis-ci.org/miguel-porto/SiMRiv)

**Individual-based simulation of multistate movements in rivers, heterogeneous and homogeneous spaces incorporating local landscape bias in R.**

Available on [CRAN](https://cran.r-project.org/web/packages/SiMRiv/index.html), can be directly installed from within R.

Provides functions to generate and analyze spatially-explicit individual-based multistate movements in rivers, heterogeneous and homogeneous spaces.
This is done by incorporating landscape bias on local behaviour, based on resistance rasters.
Although originally conceived and designed to simulate trajectories of species constrained to linear habitats/dendritic ecological networks (e.g. river networks), the simulation algorithm is built to be
highly flexible and can be applied to any (aquatic, semi-aquatic or terrestrial) organism, independently on the landscape in which it moves.
Thus, the user will be able to use the package to simulate movements either in homogeneous landscapes, heterogeneous landscapes (e.g. semi-aquatic animal moving mainly along rivers but also using
the matrix), or even in highly contrasted landscapes (e.g. fish in a river network). The algorithm and its input parameters are the same for all cases, so that results are comparable.

Simulated trajectories can then be used as mechanistic null models (Potts & Lewis 2014) to test a variety of 'Movement Ecology'
hypotheses (Nathan et al. 2008), including landscape effects (e.g. resources, infrastructures) on animal movement and species site fidelity (Powell 2000),
or for predictive purposes (e.g. road mortality risk, dispersal/connectivity).
The package should be relevant to explore a broad spectrum of ecological phenomena, such as those at the interface of animal behaviour, management, landscape and movement ecology, disease and invasive species
spread, and population dynamics.

