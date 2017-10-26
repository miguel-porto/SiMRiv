\name{speciesModel}
\alias{speciesModel}
\title{Defines a species model to adjust to a real trajectory}
\description{
  The sole purpose of this function is to be used in conjunction with \code{\link{adjustModel}}. It is used to tell the optimization
  algorithm which parameters are to be approximated, and which are constant.
}
\usage{
  speciesModel(type, perceptual.range = 0, steplength = 1
    , prob.upperbound = 0.5, max.concentration = 0.99)
}
\arguments{
  \item{type}{the type of movement to "fit". One of \code{CRW}, \code{RW.CRW}, \code{CRW.CRW}, \code{RW.CRW.sl}, \code{CRW.CRW.sl}}
  \item{perceptual.range}{the perceptual range for all states.}
  \item{steplength}{the fixed step length for fixed step length types \code{CRW}, \code{RW.CRW}, \code{CRW.CRW} or the maximum allowed value for variable step length types \code{RW.CRW.sl}
    , \code{CRW.CRW.sl}, \code{CRW.CRW.CRW.sl} and \code{CRW.RW.Rest.sl}.}
  \item{prob.upperbound}{the maximum allowed value for the state switching probabilities. The default is 0.5 because very high state switching probabilities don't make much sense from a biological point of view.}
  \item{max.concentration}{the maximum allowed value for the turning angle concentration parameter for CRWs [0, 1]. By default it is set to 0.99 because values higher than this (for technical reasons) result in straight line paths, which is a technical artifact.}
}
\details{
  This function defines the type of movement to be adjusted with \code{\link{adjustModel}}. Before choosing the type, it is good practice
  to plot the real trajectory and visually assess which would be the most adequate model to try. Currently included movement types are:

  \itemize{
  \item{\code{CRW}: single state CRW, fixed step length (1 parameter)}
  \item{\code{RW.CRW}: two state RW/CRW, fixed step length (3 parameters)}
  \item{\code{CRW.CRW}: two state CRW/CRW, fixed step length, (4 parameters)}
  \item{\code{RW.CRW.sl}: two state RW/CRW, variable step length, (5 parameters)}
  \item{\code{CRW.CRW.sl}: two state CRW/CRW, variable step length (6 parameters)}
  \item{\code{CRW.CRW.CRW.sl}: three state CRW/CRW/CRW, variable step length (12 parameters)}
  \item{\code{CRW.RW.Rest.sl}: three state CRW/RW/Rest, variable step length (7 parameters)}
  }

  However, the user can easily write any custom function for addressing other movement types, see the code for details.
}
\value{
  Returns a function that creates a species from a vector of parameter values. This function is normally used to create species
  from the \code{\link{adjustModel}} results, see examples there.
}
\seealso{
  \code{\link{adjustModel}}.
}
\examples{
library(SiMRiv)

model <- speciesModel("RW.CRW.sl")

# this shows the parameters that will be approximated
model

# this creates a species with 2 states
# RW and a CRW with correlation 0.9
# with the switching probabilities RW->CRW = 0.01, CRW->RW = 0.05
# and the step lengths RW = 15, CRW = 50.

species <- model(c(0.9, 0.01, 0.05, 15, 50))
}

