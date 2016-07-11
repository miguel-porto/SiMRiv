\name{simulate}
\alias{simulate}
\title{Simulate movements in river networks or heterogeneous landscapes}
\description{
  Fast and spatially-explicit simulation of multi-state movements of individuals created with the function \code{\link{species}} in an optional landscape resistance raster.
}
\usage{
simulate(individuals, time, coords = NULL
  , states = NULL, resist = NULL)
}
\arguments{
  \item{individuals}{a \code{species} or a list of N \code{species} whose movements are to be simulated}
  \item{time}{the number of time steps to be simulated}
  \item{coords}{a N x 2 matrix giving the initial coordinates of the simulation, for each individual. Can be a vector of the form c(x, y) if only one individual is provided.}
  \item{resist}{an optional landscape resistance raster of class
  \code{RasterLayer} (package \code{raster}). If not provided, movements are simulated in an homogeneous environment.}
  \item{states}{Not implemented yet. An optional numeric vector of length N giving the initial state of each individual.}
}
\details{
  Performs a mechanistic simulation of the movement of the given individual(s) simultaneously, in the given landscape raster defining physical resistance values.
  At present, multiple individuals do not interact, but in the upcoming version it will be possible to
  define positive and negative interactions between simulated individuals, thus accounting for spatial bias.

  The simulation runs in a series of micro-steps, and is intended to be a high-resolution simulation
  (which can be later sampled with \code{\link{sampleMovement}} to emulate real field data, e.g. telemetry data).
  
  In summary, in each micro-step the individual chooses a random direction which is based on the previous step heading and
  in the resistance context at the current position, such that it will avoid heading to areas with high resistance.
  This evaluation depends on the individual's perception window in the current movement state.
  At each step, a test to see if the individual changes its state is also performed, based on the provided \code{\link{transitionMatrix}}.
  Each state has its own perception window, step length and angular correlation (with previous step heading).

  In more detail, in each of the \code{time} steps, the procedure is as follows:
  \enumerate{
    \item{Draw state for this step according to state transition matrix and previous state}
    \item{Compute empirical probability density from the landscape raster values around current position (resistance component). See below.}
    \item{Compute probability density given the previous step heading and correlation defined in the current state (correlated walk component)}
    \item{Intersect the resistance component with the correlated walk component to make a compound probability density}
    \item{Draw the new heading from the probability density distribution computed above}
    \item{Compute the length of the step that will be taken as a fraction of current step's defined length proportional to mean resistance of the starting and ending points in the chosen heading}
    \item{Move to the new position}
  }
  \subsection{Details of the simulation algorithm}{
    \subsection{The landscape resistance raster}{
      This raster represents the amount of physical resistance that is offered to the simulated individuals. The values must be between 0 (no resistance) and 1 (infinite resistance).
      In the future, other types of rasters can be provided, for example rasters for resource
      availability, habitat suitability and points of attraction/repulsion, allowing to conduct simulations
      with various types of spatial bias.

      A careful choice of pixel size must be taken for the resistance raster. If rasterizing from vector lines (e.g. river network)
      , be sure to adjust the pixel size so that there are no gaps between river pixels and all pixels of the river are connected
      orthogonally.
    }
    \subsection{The empirical probability density}{
      The empirical probability density for a given point (resistance component) is computed by summing the (1-resistance) values along all discrete radial lines departing from that point, forming a circle.
      The length of the lines (i.e. radius) and weighting given to each pixel are defined in the current state's perception window.
      The sums are packed and used as the circular empirical distribution of the resistance component.
      This will be crossed with the correlated walk component to yield the final empirical probability distribution from which heading will be drawn.
    }
  }
}
\value{
	A matrix with 3 columns for each simulated individual, in the order x1, y1, state1, x2, y2, state2, ...; and the same number of rows as the simulation length (given by \code{time}).
}
\note{The returned object will change in the future.}
\seealso{
  \code{\link{species}}, \code{\link{sampleMovement}}.
}
\examples{
library(SiMRiv)

## A classic: simple random walk
## i.e. single-state uncorrelated movement in an homogeneous landscape
######################################################################

rand.walker <- species(state.randomWalk())   # a single state, other parameters
                                             # set to defaults
sim <- simulate(rand.walker, 10000)
plot(sim, type="l", asp=1)

## two random walkers
#####################

sim <- simulate(list(rand.walker, rand.walker), 10000)
plot(sim[,1:2], type="l", asp=1, xlim=range(sim), ylim=range(sim), col=2)
lines(sim[,4:5], col=3)

## Another classic: Levy walk-like movement
## i.e. two-state movement: random walk
## with bursts of correlated random walk
###########################################

LevyWalker <- species(
    list(state.randomWalk(),
    state.CRW(0.99)),
    transitionMatrix(0.005, 0.02))
    
sim <- simulate(LevyWalker, 10000)
plot(sim, type="l", asp=1)

## Heterogeneous landscapes: fish in a river network
####################################################

data(river)   # sample river raster in a fish's perspective, i.e.
              # resistance is 0 within the river, 1 otherwise.

## let's try a Levy-like movement in a river network

LevyWalker <- species(list(
  state(0, perceptionWindow("cir", 100), 10, "RandomWalk")
    # perception window radius (100) and step length (20) must be
    # adequate to the raster resolution!
  ,state(0.97, perceptionWindow("cir", 500), 20, "CorrelatedRW")  # idem
), transitionMatrix(0.005, 0.001))

sim <- simulate(LevyWalker, 100000, resist = river, coords = c(280635, 505236))

## plot movement; we use a high-res TIFF so that it can be viewed in detail

tiff("movement.tif", wid=10000, hei=10000, comp="lzw")
par(mar = c(0, 0, 0, 0))
plot(river, asp = 1, col = gray(seq(1, 0.5, len = 2))
  , ylim = range(sim[,2]), xlim = range(sim[,1]), axes = FALSE)
lines(sim, lwd = 1.5, col = "#0000ffcc")
dev.off()

## if we want the kernel density overlaid, uncomment these and put before dev.off()
# library(ks)
# d <- kde(sim[,1:2])
# plot(d, disp = "image", add=TRUE, col = rgb(1, 0, 0, seq(0, 1, len = 15)))
}
\keyword{simulation}

