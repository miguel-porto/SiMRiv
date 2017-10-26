\name{sampleMovement}
\alias{sampleMovement}
\title{Resample a simulated movement and compute step-wise statistics}
\description{
  Resamples a movement simulated with the function \code{\link{simulate}} to a lower temporal resolution (frequency).
  Simultaneously, computes step length, turning angles (\cite{Turchin 1998}) and accumulated resistance
  for each resampled step, assuming that a straight line connects each re-sampled location.
}
\usage{
  sampleMovement(relocs, resolution = 1, resist = NULL)
}
\arguments{
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
