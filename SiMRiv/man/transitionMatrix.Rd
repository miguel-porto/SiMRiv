\name{transitionMatrix}
\alias{transitionMatrix}
\title{Define a state transition matrix}
\description{
  Defines a state transition matrix to be used when creating \code{\link{species}}.
}
\usage{
  transitionMatrix(...)
}
\arguments{
  \item{...}{the probabilities that will form the matrix, see Details. If none given, returns a 1-element matrix (for one state only)}
}
\details{
  The transition matrix (Markov matrix) is a square, non-symmetric matrix with all elements between 0 and 1, and whose rows must sum to 1 (but not columns).
  It defines the probability of the individual changing from each behavioral state to another, and this
  is tested in each time step of the simulation, hence probabilities should be small.
  
  This function is just a helper to create such matrix.
  The arguments are probabilities given in the following order (example for 3 states):

  Probability of changing from:

  state 1 --> state 2
  
  state 1 --> state 3
  
  state 2 --> state 1
  
  state 2 --> state 3
  
  state 3 --> state 1
  
  state 3 --> state 2


  The diagonal (probablity of remaining in the same state) is computed so that rows sum to 1.  
}
\value{
  A numeric matrix.
}
\seealso{
  \code{\link{species}}.
}
\examples{
## a 3-state transition matrix

transitionMatrix(0.01,0.02,0,0.03,0.0001,0)
}
