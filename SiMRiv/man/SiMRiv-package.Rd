\name{SiMRiv-package}
\alias{SiMRiv-package}
\alias{SiMRiv}
\docType{package}
\title{
\packageTitle{SiMRiv}
}
\description{
Provides functions to generate and analyze individual-based, spatially-explicit simulations of multi-state movements in heterogeneous
landscapes, based on "resistance" rasters. Although originally conceived and designed to simulate spatially-explicit
trajectories of species constrained to linear, dendritic habitats (e.g., river networks), the simulation algorithm is built to be highly
flexible and can be applied to any (aquatic, semi-aquatic or terrestrial) organism. Thus, the user will be able to use the package to
simulate movements either in homogeneous landscapes, heterogeneous landscapes (e.g. semi-aquatic animal moving mainly along rivers but also using the matrix),
or even in highly contrasted landscapes (e.g. fish in a river network). The algorithm and its input parameters are the same for all cases, so that
results are comparable. Simulated trajectories can then be used as mechanistic null models (Moorcroft and Lewis 2006) to test e.g. for
species site fidelity (Powell 2000) and other 'Movement Ecology' hypotheses (Nathan et al. 2008), or for other predictive purposes.
The package should thus be relevant to explore a broad spectrum of ecological phenomena, such as those at the interface of animal behaviour,
landscape, spatial and movement ecology, disease and invasive species spread, and population dynamics.
}
\details{
\code{\link{simulate}} is the central function. See the examples in \code{?simulate} to quickly get started, or the vignette for a more verbose tutorial.

\packageIndices{SiMRiv}
}
\references{
  \itemize{
    \item{Powell, R. A. 2000. Animal home ranges and territories and home range estimators. In: Research techniques in animal ecology: controversies and consequences, 442. Boitani, L., & Fuller, T. (Eds.). Columbia university press, New York: pp.65-110.}
    \item{Moorcroft, P. R. & Lewis, M. A. 2006. Mechanistic Home Range Analysis. Monographs in Population Biology 43. Eds. Levin S.A. and H.S. Horn. Princeton University Press. pp 172.}
    \item{Nathan, R., Getz, W. M., Revilla, E., Holyoak, M., Kadmon, R., Saltz, D., & Smouse, P. E. 2008. A movement ecology paradigm for unifying organismal movement research. Proceedings of the National Academy of Sciences,105(49), 19052-19059.}
  }
}

\keyword{ package }

\examples{
## a simple Levy-like movement in homogeneous space
## see ?simulate for more complex examples

LevyWalker <- species(
    state.RW() + state.CRW(0.99),
    transitionMatrix(0.005, 0.02))

sim <- simulate(LevyWalker, 20000)
plot(sim, type="l", asp=1)
}
