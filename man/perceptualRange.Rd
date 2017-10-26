\name{perceptualRange}
\alias{perceptualRange}
\title{Define a perceptual range}
\description{
  Defines the perceptual range to be used in a movement state.
}
\usage{
perceptualRange(type = "circular", radius)
}
\arguments{
  \item{type}{defines the type of weights that are given to each pixel, according to the distance to its center. One of \code{circular} or \code{gaussian}.}
  \item{radius}{the radius of the circular perceptual range, or the sigma of the gaussian perceptual range, in map units.}
}
\details{
  The perceptual range is often defined as the distance (radius) at which the species perceives the environment, 
  based, e.g., on sense of smell, vision, audition, etc. (\cite{Lima & Zollner 1996}; \cite{Powell 2000}).
  In SiMRiv, perceptual range should be seen as the distance (radius) at which the species evaluates the landscape resistance around its current location, influencing species next heading.
  Perceptual range size can be defined (in meters or other map units) based on available literature (on species perceptual range, or, as its surrogate, species home range size), on expert-based criteria, or be estimated from real data.
  
  A \code{circular} range gives equal weight to all pixels, which form a circle centered on current individual's position.
  A \code{gaussian} range gives weights corresponding to a gaussian kernel centered on current individual's position.
}
\seealso{
  \code{\link{state}}.
}
\references{
  \itemize{
   \item{Lima, S. L., & Zollner, P. A. (1996). Towards a behavioral ecology of ecological landscapes. Trends in Ecology & Evolution, 11(3), 131-135.}
   \item{Powell, R. A. 2000. Animal home ranges and territories and home range estimators. In: Research techniques in animal ecology: controversies and consequences, 442. Boitani, L., & Fuller, T. (Eds.). Columbia university press, New York: pp.65-110.}
  }
}
