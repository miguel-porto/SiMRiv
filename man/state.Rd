\name{state}
\alias{state}
\alias{state.Resting}
\alias{state.RW}
\alias{state.CRW}
\title{Define a movement state}
\description{
  Defines a behavioral state to be used when creating \code{\link{species}}.
}
\usage{
  state(concentration, pwind = perceptualRange("circular", 0)
    , steplen = 1, name = "")
  state.Resting()    # still state
  state.RW()         # uniform random walk (brownian motion),
                     # independent of resistance
  state.CRW(concentration)  # correlated random walk,
                          # independent of resistance
}
\arguments{
  \item{concentration}{turning angle concentration, a value between 0 (uniform distribution resulting in random walk) and 1 (only one value possible resulting in a straight line path)}
  \item{pwind}{a \code{\link{perceptualRange}} definition}
  \item{steplen}{the base (maximum) step length of this state in map units. Note that the actual step length depends on the resistance in each step.}
  \item{name}{the name of the state}
}
\details{
  See \code{\link{Arith-methods}} for more convenient ways of setting parameters, instead of using \code{state}.
}
\value{
  An object of class \code{state}.
}
\note{
  The perceptual range radius and step length must be adequate to the resolution of the resistance raster (if provided in simulations).
  If no raster will be provided, then the perceptual range is irrelevant, and the step length has solely a relative meaning (in relation to other states or other species).
  
  For a review of different random walks, see Codling et al. (2008)
}
\seealso{
  \code{\link{species}}, \code{\link{perceptualRange}}, \code{\link{Arith-methods}}.
}
\examples{
## a correlated random walk influenced by landscape

state(0.97, perceptualRange("cir", 500), 10, "CorrelatedRW")

## the same, but using the shortcut form

state.CRW(0.97) * 500 + 10
}
\references{
  \itemize{
    \item{Codling, E. A., Plank, M. J., & Benhamou, S. 2008. Random walk models in biology. Journal of the Royal Society Interface, 5(25), 813-834.}
  }
}
