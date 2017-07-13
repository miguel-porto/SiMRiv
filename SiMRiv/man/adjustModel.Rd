\name{adjustModel}
\alias{adjustModel}
\title{Finds ("estimates") simulation input parameters able to replicate a given (real) trajectory, assuming the given species model}
\description{
  Given a trajectory, a type of movement and the time resolution at which the user wants to simulate, 
  this function approximates the values for the input parameters so that the simulated movement is maximally
  similar to the given trajectory, in terms of general non-spatial patterns. If the user wants to simulate at a higher frequency 
  than real data (which is the norm), the function tries to maximize the similarity between the real trajectory and 
  the simulated trajectories after downsampled to the same frequency as the real.
}
\usage{
  adjustModel(realData, species.model, resolution = 10
    , resistance = NULL, coords = NULL, start.resistance = NULL
    , nbins = 100, window.size = dim(reference$stats)[1] %/% nbins
    , nbins.hist = c(5, 5), step.hist.range = c(0, 1)
    , step.hist.log = FALSE, nrepetitions = 1
    , popsize = 100, generations = seq(5, 1000, by=5), mprob = 0.2
    , parallel = is.null(resistance)
  )
}
\arguments{
  \item{realData}{the given trajectory for which the simulation input parameters are to be "estimated", given as a matrix with two columns (coordinates) and assuming that relocations are equally spaced in time.}
  \item{species.model}{the species model to adjust, created with \code{\link{speciesModel}}. This defines the type of movement that is to be adjusted (e.g. how many, and which type of, behavioral states, see details).}
  \item{resolution}{the desired time frequency of the simulations for which parameters will be approximated, as a fraction of the real data, i.e. a value of 20 will simulate movements at a 20-fold higher frequency than real data.}
  \item{resistance}{the resistance raster to use in simulations during parameter approximation.}
  \item{coords}{the initial coordinates of the simulated individuals (only relevant if \code{resistance} is provided because the metrics used in optimization are spatially-agnostic).}
  \item{start.resistance}{the maximum value of resistance in which the individuals are allowed to start. Only used if \code{resistance} is provided and \code{coords} is not provided.}
  \item{nbins}{the number of time lags within which the standard deviation of the turning angles will be computed during optimization. Ignored if \code{window.size} is provided.}
  \item{window.size}{the size (in steps) of the time lags within which the standard deviation of the turning angles will be computed during optimization. A different way of providing \code{nbins}. See details.}
  \item{nbins.hist}{a vector with two positive integers defining the number of histogram bins for turning angle variation and step length histograms respectively. These bins will be used during optimization to compare simulated trajectories with the real trajectory to infer the quality of the "fit". Comparison is done bin-wise, which means that the more bins, the more detailed comparison is, the higher the quality of the solutions, but the more difficult it is to find them.}
  \item{step.hist.range}{the quantiles used to define the range of the step length histogram computation. Used if the user wants to exclude outliers. The default is [min, max].}
  \item{step.hist.log}{set to \code{TRUE} to use the histogram of the logarithm of the step lengths rather than of the raw values in the comparisons. Setting to \code{TRUE} usually results in more detail in the comparisons.}
  \item{nrepetitions}{the number of simulations conducted for each solution evaluation during optimization. If >1, the quality of the solutions is computed by comparing the averaged histograms across repetitions with the real histograms.}
  \item{popsize}{optimization algorithm option, to pass to \code{\link{nsga2}}}
  \item{generations}{optimization algorithm option, to pass to \code{\link{nsga2}}}
  \item{mprob}{optimization algorithm option, to pass to \code{\link{nsga2}}}
  \item{parallel}{set to \code{TRUE} to use multicore processing.}
}
\details{
  This function finds possible parameters for the simulation (solutions), so that the resulting movements are as similar as possible
  to the given real trajectory, in terms of their intrinsic properties measured by step lengths and variation in turning angles. The input parameter
  approximations are found using a multiobjective genetic algorithm (NSGA-II, \cite{Deb et al. 2002}). The algorithm minimizes a vector
  of objectives (\code{sum(nbins.hist)}) whose values are computed by the absolute differences between each pair of bins (real and simulated)
  of two histograms: an histogram of the step lengths and an histogram of the standard deviation in turning angles computed
  in a moving time window along each trajectory.

  SiMRiv simulations are intended to reproduce the fine-scale movement steps, unlike what is normally collected in field data.
  Hence, simulations should be conducted with a much higher time frequency than provided by the real data, which poses challenges for parameterization.
  This function incorporates this difference in the time scale during optimization (\code{resolution}), allowing the user
  to find the input parameters for simulations at a much higher frequency, which, when downsampled to the real data's time
  frequency, will present similar patterns.
  
  There are no limits to the number of parameters that can be approximated, but obviously, the higher the number, the more difficult convergence
  is achieved. As in any other method, a compromise should be sought. The number of parameters to approximate is defined by the user by providing a
  \code{\link{speciesModel}}. This defines how many states and which types of states are to be "fit" to the data. See \code{\link{speciesModel}} for details.
  A good starting point is to provide a two-state species model, in which both states are Correlated Random Walks. This model involves the approximation
  of 6 parameters and is sufficiently flexible for simulating a variety of movements, while not overly complex. Note that all complex models can
  accomodate to simpler ones (i.e. the simpler models are special cases of the complex ones).
}
\note{
  This function is an experimental feature, here provided only to guide the user on how to parameterize the simulations. Care must be taken
  when interpreting the results, specially by assessing algorithm convergence and visually comparing simulations with the approximated
  parameters to the real data.
}
\value{
  The object returned by \code{\link{nsga2}} (package \code{mco}), see details therein. If \code{generations} is a vector (which is recommended,
  for assessing convergence), this object contains the approximated input parameter values in each generation given in the vector. See examples
  for easily plotting results.
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
