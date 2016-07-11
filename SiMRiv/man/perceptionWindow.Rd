\name{perceptionWindow}
\alias{perceptionWindow}
\title{Define a perception window}
\description{
  Defines the perception window to be used in a movement state.
}
\usage{
perceptionWindow(type = "circular", radius)
}
\arguments{
  \item{type}{defines the type of weights that are given to each pixel, according to the distance to its center. One of \code{circular} or \code{gaussian}.}
  \item{radius}{the radius of the circular perception window, or the sigma of the gaussian perception window, in map units.}
}
\details{
  The perception window is defined as the distance radius at which the species perceives the environment,
  based e.g. on sense of smell, vision, audition, etc. (\cite{Powell 2000}).
  
  A \code{circular} window gives equal weight to all pixels, which form a circle centered on current individual's position.
  
  A \code{gaussian} window gives weights corresponding to a gaussian kernel centered on current individual's position.
}
\seealso{
  \code{\link{state}}.
}
\references{
  \itemize{
    \item{Powell, R. A. 2000. Animal home ranges and territories and home range estimators. In: Research techniques in animal ecology: controversies and consequences, 442. Boitani, L., & Fuller, T. (Eds.). Columbia university press, New York: pp.65-110.}
  }
}
