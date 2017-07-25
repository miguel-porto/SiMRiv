\name{species}
\alias{species}
\title{Create a species}
\description{
  Creates a species, characterized by one or more behavioral states, to be simulated with function \code{\link{simulate}}.
}
\usage{
  species(states, trans = transitionMatrix(), name = "<unnamed>", resistanceMap = NULL)
}
\arguments{
  \item{states}{a list of \code{state}s characterizing the behavior of the species, or a single \code{state}, for simple movements}
  \item{trans}{a square state transition matrix, defining the probability of changing between states. For convenience, use the function \code{\link{transitionMatrix}}. Can be omitted if this species has a single-state movement.}
  \item{name}{the name of the species}
  \item{resistanceMap}{not used. Will be implemented in future versions.}
}
\details{
  The rows and columns of the transition matrix correspond in the same order to the list of states.
  The matrix is not symmetric, and is read along the rows, i.e. the probability of changing
  from state 2 to state 1 is located in row 2, column 1; hence rows must sum to 1 but columns not.
}
\value{
  An object of class \code{species}.
}
\seealso{
  \code{\link{simulate}}, \code{\link{perceptualRange}}, \code{\link{state}}, \code{\link{transitionMatrix}}, \code{\link{Arith-methods}}.
}
\examples{
## example from 'simulate'

## note: perceptual range radii and step lengths must be
## adequate to the raster resolution!

LevyWalker <- species(
	(state.RW() * 100 + 10) + (state.CRW(0.97) * 500 + 20)
, transitionMatrix(0.005, 0.001))
}
