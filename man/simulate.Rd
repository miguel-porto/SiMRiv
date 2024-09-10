\name{simulate}
\alias{simulate}
\title{Simulate movements in river networks, homogeneous, or heterogeneous landscapes}
\description{
  Performs fast and spatially-explicit simulation of multi-state random movements (\cite{Morales et al. 2004, McClintock et al. 2012}) of individuals created with the function \code{\link{species}} in an optional landscape resistance raster.
}
\usage{
simulate(individuals, time, coords = NULL
  , states = NULL, resist = NULL, angles = NULL
  , start.resistance)
}
\arguments{
  \item{individuals}{a \code{species} or a list of N \code{species} whose movements are to be simulated}
  \item{time}{the number of time steps to be simulated}
  \item{coords}{a N x 2 matrix giving the initial coordinates of the simulation, for each individual. Can be a vector of the form c(x, y) if only one individual is provided.}
  \item{resist}{an optional landscape resistance raster of class \code{RasterLayer}. If not provided, movements are simulated in an homogeneous environment.}
  \item{angles}{an optional numeric vector of length N defining the initial heading of each individual, in radians. Zero is north and angles increase clockwise.}
  \item{start.resistance}{an optional scalar in the range [0, 1] giving the maximum resistance value in which the individuals are allowed to start, if \code{coords} is \code{NULL} and \code{resist} is provided}
  \item{states}{Not implemented yet. An optional numeric vector of length N defining the initial state of each individual.}
}
\details{
  Performs a mechanistic simulation of the movement of the given individual/s (when more than one individual is given,
  their movements are simulated simultaneously) in the given landscape raster defining physical resistance values.
  At present, multiple individuals do not interact, but in the upcoming version it will be possible to
  define positive and negative interactions between simulated individuals, thus accounting for spatial bias.

  The simulation runs in a series of micro-steps, and is intended to be a high-resolution simulation
  (which can be later sampled with \code{\link{sampleMovement}} to emulate real field data, e.g. telemetry data).
  
  In summary, at each micro-step the individual chooses a random direction which is based on the previous step heading and
  in the resistance context at the current position, such that it will avoid heading to areas with high resistance.
  This evaluation depends on the individual's perceptual range in the current movement state.
  At each step, a test to see if the individual changes its state is also performed, based on the provided \code{\link{transitionMatrix}}.
  Each state may have its own perceptual range, step length and angular correlation (with previous step heading). It's up to the user the
  definition of these values (e.g., expert- or literature-based), but we provide an experimental function to numerically approximate these
  values from real data, see \code{\link{adjustModel}}.

  In more detail, in each of the \code{time} steps, the procedure is as follows:
  \enumerate{
    \item{Draw state for the current step according to state transition matrix and previous state}
    \item{Compute empirical probability density for changing heading, from the landscape raster values around current position (resistance component). See details below.}
    \item{Compute probability density for changing heading, given the previous step heading and correlation defined in the current state (correlated walk component)}
    \item{Intersect the resistance component with the correlated walk component to make a compound probability density for changing heading}
    \item{Draw the new heading from the probability density distribution computed above}
    \item{Compute the length of the step that will be taken as a fraction of current step's defined length proportional to mean resistance of the starting and ending points in the chosen heading}
    \item{Move to the new position, defined by the drawn heading and length of step.}
  }
  \subsection{Details of the simulation algorithm}{
    \subsection{The landscape resistance raster}{
      This raster represents the amount of physical resistance that is offered to the simulated individuals. The values must be between 0 (no resistance) and 1 (infinite resistance).
      In the future, other types of rasters can be provided, for example rasters for resource
      availability, habitat suitability and points of attraction/repulsion, allowing to conduct simulations
      with various types of spatial bias.

      A careful choice of pixel size must be taken for the resistance raster. If rasterizing from vector lines (e.g. river network),
      please be sure to adjust the pixel size so that there are no gaps between river pixels and all pixels of the river are connected
      orthogonally.
    }
    \subsection{The empirical probability density}{
      The empirical probability density for a given point (resistance component) is computed by summing the \code{1 - resistance} values along a set of discrete radial lines departing from that point, forming a circle.
      The length of the lines (i.e. radius) and weighting given to each pixel are defined in the current state's perceptual range.
      The sums are packed and used as the circular empirical distribution of the resistance component.
      This will be crossed with the correlated walk component to yield the final empirical probability distribution from which heading will be drawn.
    }
  }
}
\value{
	A matrix with 3 columns for each simulated individual, in the order x1, y1, state1, x2, y2, state2, ...; and the same number of rows as the simulation length (given by \code{time}).
}
\note{The structure of the returned object will change in the upcoming version.}
\seealso{
  \code{\link{species}}, \code{\link{sampleMovement}}, \code{\link{adjustModel}}.
}
\examples{
library(SiMRiv)

## A classic: simple random walk (Brownian motion) (Turchin 1998)
## i.e. single-state uncorrelated movement in an homogeneous landscape
######################################################################

## a single state, other parameters set to defaults

rand.walker <- species(state.RW())
sim <- simulate(rand.walker, 10000)
plot(sim, type="l", asp=1)

## two random walkers
#####################

sim <- simulate(list(rand.walker, rand.walker), 10000)
plot(sim[,1:2], type="l", asp=1, xlim=range(sim), ylim=range(sim), col=2)
lines(sim[,4:5], col=3)

## Another classic: Levy walk-like movement (e.g. Sims et al. 2008)
## i.e. two-state movement: composition of small-scale random walks
## with bursts of longer, correlated random walks
###################################################################

LevyWalker <- species(
    state.RW() + state.CRW(0.99),
    transitionMatrix(0.005, 0.02))
    
sim <- simulate(LevyWalker, 10000)
plot(sim, type="l", asp=1)

## Linear habitats, e.g. fish in a river network
################################################

## load sample river raster in a fish's perspective,
## i.e. resistance is 0 within the river, 1 otherwise.

river <- terra::rast(system.file("doc/river.tif", package="SiMRiv"))

## let's try a Levy-like movement in a river network
## note: perceptual range radii and step lengths must be
## adequate to the raster resolution!

LevyWalker <- species(list(
    state(0, perceptualRange("cir", 100), 10, "RandomWalk")
    ,state(0.97, perceptualRange("cir", 500), 20, "CorrelatedRW")
), transitionMatrix(0.005, 0.001))

## NOTE: the following lines do exactly the same as above, but
## using the more convenient arithmetic operator shortcuts

LevyWalker <- species(
	(state.RW() * 100 + 10) + (state.CRW(0.97) * 500 + 20)
, transitionMatrix(0.005, 0.001))

sim <- simulate(LevyWalker, 20000, resist = river
    , coords = c(280635, 505236))

## plot movement; we use a high-res TIFF so that it
## can be viewed in detail
\dontrun{
tiff("movement.tif", wid=5000, hei=5000, comp="lzw")
par(mar = c(0, 0, 0, 0))
plot(river, asp = 1, col = gray(seq(1, 0.5, len = 2))
    , ylim = range(sim[,2]), xlim = range(sim[,1]), axes = FALSE)
lines(sim, lwd = 2, col = "#0000ffcc")
dev.off()
}

## if we want the kernel density overlaid,
## uncomment these and put before dev.off()
# library(ks)
# d <- kde(sim[,1:2])
# plot(d, disp = "image", add=TRUE
#    , col = rgb(1, 0, 0, seq(0, 1, len = 15)))
}
\references{
  \itemize{
    \item{McClintock, B. T., King, R., Thomas, L., Matthiopoulos, J., McConnell, B. J., & Morales, J. M. 2012. A general discrete-time modeling framework for animal movement using multistate random walks. Ecological Monographs, 82(3), 335-349.}
    \item{Morales, J. M., Haydon, D. T., Frair, J., Holsinger, K. E., & Fryxell, J. M. 2004. Extracting more out of relocation data: building movement models as mixtures of random walks. Ecology, 85(9), 2436-2445.}
    \item{Sims, D. W., Southall, E. J., Humphries, N. E., Hays, G. C., Bradshaw, C. J., Pitchford, J. W., ... & Morritt, D. 2008. Scaling laws of marine predator search behaviour. Nature, 451(7182), 1098-1102.}
    \item{Turchin, P. 1998. Quantitative analysis of movement: measuring and modeling population redistribution in animals and plants (Vol. 1). Sinauer Associates, Sunderland, MA.}
  }
}
\keyword{simulation}

