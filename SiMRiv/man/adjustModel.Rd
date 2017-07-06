\name{adjustModel}
\alias{adjustModel}
\title{Finds ("estimates") simulation input parameters able to replicate a given (real) trajectory, assuming the given species model}
\description{
  Given a trajectory, a type of movement and the time resolution at which the user wants to simulate, 
  this function approximates the values for the input parameters so that the simulated movement is maximally
  similar to the given trajectory, in terms of general non-spatial patterns. If the user wants to simulate at a higher frequency 
  than real data (which is the norm), the function tries to maximize the similarity between the real trajectory and 
  the simulated trajectories after downsampled to the same frequency as the real.
}
\usage{
  adjustModel(realData, species.model, resolution = 10
    , resistance = NULL, coords = NULL, start.resistance = NULL
    , popsize = 100, ngenerations = 400, mprob = 0.2
    , nbins = 100, nbins.hist = 7, nrepetitions = 1
    , parallel = is.null(resistance)
  )
}
\arguments{
  \item{realData}{the given trajectory for which the simulation input parameters are to be "estimated", given as a matrix with two columns (coordinates) and assuming that relocations are equally spaced in time.}
  \item{species.model}{the species model to adjust, created with \code{\link{speciesModel}}. This defines the type of movement that is to be adjusted (e.g. how many, and which type of, behavioral states, see details).}
  \item{resolution}{the desired time frequency of the simulations for which parameters will be approximated, as a fraction of the real data, i.e. a value of 20 will simulate movements at a 20-fold higher frequency than real data.}
  \item{resistance}{the resistance raster to use in simulations during parameter approximation.}
  \item{coords}{the initial coordinates of the simulated individuals (only relevant if \code{resistance} is provided because the metrics used in optimization are spatially-agnostic).}
  \item{start.resistance}{the maximum value of resistance in which the individuals are allowed to start. Only used if \code{resistance} is provided and \code{coords} is not provided.}
  \item{nbins}{the time lag used for computing the standard deviation of the turning angles}
	, start.resistance = NULL
# fitness function parameters
	, nbins = 100
	, nbins.hist = 7
	, nrepetitions = 1
# GA options
	, popsize = 100, ngenerations = 400, mprob = 0.2
	, parallel = is.null(resistance)

  \item{relocs}{the simulated movement, an object returned by \code{\link{simulate}}. NOTE: currently this is only implemented for simulations of single individuals.}
  \item{resolution}{movement will be resampled every this number of time ticks. If 1, no resampling is done (but metrics are computed).}
  \item{resist}{a landscape resistance raster, usually the same that was used in \code{simulate}}
}
\details{
  This function mimics what happens in real world movement data: it resamples the simulated movement
  (which is supposed to be infinitesimal) into a lower temporal resolution, so that it is comparable to real world field data (e.g. telemetry data).
  During the process, it computes, for each resampled step:
  \itemize{
  \item{the step length}
  \item{the turning angle}
  \item{the accumulated resistance along the step, assuming a straight line is taken from start to end}
  }
}
\note{
  These metrics are only meaningful for \code{resolution} >> 1, otherwise they are just a consequence of the simulation input parameters.
}
\value{
  A list with the components \code{relocs} and \code{stats}.
  
  \code{relocs} contains the resampled positions, \code{stats} contains the metrics for each step (which has N-2 rows because of the turning angles).
}
\seealso{
  \code{\link{simulate}}.
}
\examples{
library(SiMRiv)

LevyWalker <- species(
  state.RW() + state.CRW(0.99),
  trans = transitionMatrix(0.005, 0.02))

sim <- simulate(LevyWalker, 10000)
resamp <- sampleMovement(sim, 50)
plot(sim, type="l", asp=1, col = "#777777")
lines(resamp$relocs, col = "red")
}
\references{
  \itemize{
    \item{Turchin, P. 1998. Quantitative analysis of movement: measuring and modeling population redistribution in animals and plants (Vol. 1). Sinauer Associates, Sunderland, MA.}
  }
}
