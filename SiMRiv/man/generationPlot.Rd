\name{generationPlot}
\alias{generationPlot}
\title{Finds ("estimates") simulation input parameters able to replicate a given (real) trajectory, assuming the given species model}
\description{
  Given a trajectory, a type of movement and the time resolution at which the user wants to simulate, 
  this function approximates the values for the input parameters so that the simulated movement is maximally
  similar to the given trajectory, in terms of general non-spatial patterns. If the user wants to simulate at a higher frequency 
  than real data (which is the norm), the function tries to maximize the similarity between the real trajectory and 
  the simulated trajectories (at a higher frequency) after downsampled to the same frequency as the real.
}
\usage{
  adjustModel(realData, species.model, resolution = 10
    , resistance = NULL, coords = NULL, start.resistance = NULL
    , nbins = 100, window.size = dim(reference$stats)[1]\%/\%nbins
    , nbins.hist = c(3, 3), step.hist.range = c(0, 1)
    , step.hist.log = FALSE, nrepetitions = 1
    , popsize = 100, generations = seq(5, 1000, by=5), mprob = 0.2
    , parallel = is.null(resistance), trace = TRUE
  )
}
\arguments{
  \item{realData}{the given trajectory for which the simulation input parameters are to be "estimated", given as a matrix with two columns (coordinates) and assuming that relocations are equally spaced in time.}
  \item{species.model}{the species model to adjust, created with \code{\link{speciesModel}}. This defines the type of movement that is to be adjusted (e.g. how many, and which type of, behavioral states, see details).}
  \item{resolution}{the desired time frequency of the simulations for which parameters will be approximated, as a fraction of the real data, i.e. a value of 20 will simulate movements at a 20-fold higher frequency than real data.}
  \item{resistance}{the resistance raster to use in simulations during parameter approximation.}
}
\details{
  The generation plot is a good way to assess whether the final solutions converged.
  However, note that depending on how the problem is formulated, the solutions \strong{might} not be expected to converge to a single solution "type":
  for example, there can be two types of solutions, depicting two different ways of achieving the same result.
}
\note{
}
\value{
}
\seealso{
  \code{\link{simulate}}.
}
\examples{
}
\references{
}
