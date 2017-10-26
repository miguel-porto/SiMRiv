\name{Arith-methods}
\docType{methods}
\alias{Arith-methods}
\alias{+-methods}
\alias{+,list,state-method}
\alias{+,species,ANY-method}
\alias{+,state,list-method}
\alias{+,state,numeric-method}
\alias{+,state,state-method}
\alias{*,state,numeric-method}
\alias{*,species,ANY-method}
\title{Shortcuts for defining species movement states}

\description{
Shortcuts for defining species movement states with the arithmetic operators \code{ +, *}.

The \code{+} applied to states or species and a number defines the step length of one behavioural state or of all states of a species.

The \code{*} applied to states or species and a number defines the radius of the perceptual range of one behavioural state or of all states of a species.

Further, states can be combined with \code{+} to define multistate movements in a species. See examples.
}

\section{Methods}{
\describe{

\item{\code{signature(e1 = "list", e2 = "state")}}{
Adds one state to a list of states.
}

\item{\code{signature(e1 = "species", e2 = "ANY")}}{
Applies the operation to all states in the given species.
}

\item{\code{signature(e1 = "state", e2 = "list")}}{
Adds one state to a list of states.
}

\item{\code{signature(e1 = "state", e2 = "numeric")}}{
Sets the step length of the given state to the given value.
}

\item{\code{signature(e1 = "state", e2 = "state")}}{
Combines two states in a list of states to be used with \code{\link{species}}.
}
}}

\seealso{
  \code{\link{species}}, \code{\link{state}}.
}

\examples{
# define a species with two states
twostatespecies <- species(state.RW() + state.CRW(0.9)
  , transitionMatrix(0.01, 0.02))

# set the perception window radius of both states to 200
# and the step length to 10
twostatespecies <- twostatespecies * 200 + 10
}

\keyword{methods}
\keyword{math}
\keyword{spatial}
