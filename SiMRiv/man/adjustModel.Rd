\name{adjustModel}
\alias{adjustModel}
\title{Finds ("estimates") simulation input parameters able to replicate a given (real) trajectory, assuming the given species model}
\description{
  Given a trajectory, a type of movement and the time resolution at which the user wants to simulate, 
  this function approximates the values for the simulation input parameters so that the simulated movement is maximally
  similar to the given trajectory, in terms of general non-spatial patterns. If the user wants to simulate at a higher frequency 
  than real data (which is the norm), the function maximizes the similarity between the real trajectory and 
  the simulated trajectories (at a higher frequency) after downsampled to the same frequency as the real.
  It does so by running a genetic optimization algorithm.
}
\usage{
  adjustModel(realData, species.model, resolution = 10
    , resistance = NULL, coords = NULL, start.resistance = NULL
    , nbins = 100, window.size = dim(reference$stats)[1]\%/\%nbins
    , nbins.hist = c(3, 3), step.hist.range = c(0, 1)
    , step.hist.log = FALSE, nrepetitions = 1
    , popsize = 100, generations = seq(5, 1000, by=5), mprob = 0.2
    , parallel = is.null(resistance), trace = TRUE
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
  \item{window.size}{the size (in steps of the real sampling frequency) of the time lags within which the standard deviation of the turning angles will be computed during optimization. A different way of providing \code{nbins}. See details.}
  \item{nbins.hist}{a vector with two positive integers defining the number of histogram bins for turning angle variation and step length histograms respectively. These bins will be used during optimization to compare simulated trajectories with the real trajectory to infer the quality of the "fit". Comparison is done bin-wise, which means that the more bins, the more detailed comparison is, the higher the quality of the solutions, but the more difficult it is to find them. The default \code{c(3, 3)} yields good results in general.}
  \item{step.hist.range}{the quantiles used to define the range of the step length histogram computation. Used if the user wants to exclude outliers. The default is not to exclude outliers, thus \code{c(0, 1)}.}
  \item{step.hist.log}{set to \code{TRUE} to use the histogram of the logarithm of the step lengths rather than of the raw values in the comparisons. Setting to \code{TRUE} usually results in more detail in the comparisons.}
  \item{nrepetitions}{the number of simulations conducted for each solution evaluation during optimization. If >1, the quality of the solutions is computed by comparing the averaged histograms across repetitions, with the real histograms.}
  \item{popsize}{number of solutions to optimize, to pass to \code{\link{nsga2}}}
  \item{generations}{number of algorithm generations to run, to pass to \code{\link{nsga2}}. The default is a vector, so that convergence of results can be assessed along generations.}
  \item{mprob}{mutation probability, to pass to \code{\link{nsga2}}}
  \item{parallel}{set to \code{TRUE} to use multicore processing.}
  \item{trace}{set to \code{TRUE} to print the matrix of optimization objectives (rows) for each solution (columns) in each generation along optimization. These are the values that are being internally minimized. The number of rows is \code{sum(nbins.hist)}; the first \code{nbins.hist[1]} are the turning angle variation objectives, the last are the step length objectives.}
}
\details{
  This function finds possible parameters for the simulation (solutions), so that the resulting movements are as similar as possible
  to the given real trajectory, in terms of their intrinsic properties measured by step lengths and variation in turning angles. The input parameter
  approximations are found using a multiobjective genetic algorithm (NSGA-II, \cite{Deb et al. 2002}). The algorithm minimizes a vector
  of N objectives (\code{N=sum(nbins.hist)}) whose values are computed by the absolute differences between each pair of bins (real and simulated)
  of two histograms: an histogram of the step lengths and an histogram of the standard deviation in turning angles computed
  in a moving time window along each trajectory.

  SiMRiv simulations are intended to reproduce the fine-scale movement steps, unlike what is normally collected in field data.
  Hence, simulations should be conducted with a much higher time frequency than provided by the real data, which poses challenges for parameterization.
  This function incorporates this difference in the time scale during optimization (\code{resolution}), allowing the user
  to find the input parameters for simulations at a much higher frequency, which, when downsampled to the real data's time
  frequency, will present similar patterns. The higher the frequency, the more flexibility the model has to adjust to real data,
  but the more difficult it is to find solutions due to the higher impact of stochasticity.
  
  There are no limits to the number of parameters that can be approximated, but obviously, the higher the number, the more difficult convergence
  is achieved. As in any other method, a compromise should be sought. The number of parameters to approximate is defined by the user by providing a
  \code{\link{speciesModel}}. This defines how many states and which types of states are to be "fit" to the data. See \code{\link{speciesModel}} for details.
  A good starting point is to provide a two-state species model, in which both states are Correlated Random Walks. This model involves the approximation
  of 6 parameters and is sufficiently flexible for simulating a variety of movements, while not overly complex. Note that all complex models can
  accomodate to simpler ones (i.e. the simpler models are special cases of the complex ones).
}
\note{
  This function is an experimental feature, here provided only to guide the user on how to parameterize the simulations. Care must be taken
  when interpreting the results, at least by assessing algorithm convergence and visually comparing simulations with the approximated
  parameters to the real data (see examples).
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
\dontrun{
library(SiMRiv)
library(adehabitatLT)

## simulate "real" data, a Levy walk, for which we want to
## parameterize our model

real.data <- simm.levy(1:500)[[1]][, 1:2]

## Define a species model to adjust. Let's assume we don't know
## much about what kind of real data we have, hence define
## a flexible model: a two-state correlated random walk model
## with variable step lengths.
## This model implies "estimating" 6 parameters:
##   - turning angle correlation of state 1 [0, 1]
##   - turning angle correlation of state 2 [0, 1]
##   - switching probability S1 -> S2 [0, 1]
##   - switching probability S2 -> S1 [0, 1]
##   - maximum step length of state 1 [0, ?]
##   - maximum step length of state 2 [0, ?]

## Let's assume we want to simulate at a 20 times higher time frequency
## than real data.
## In order to allow our model to adjust to real data, we have
## to provide a maximum allowable step length to the optimization algorithm
## that allows to recover real data after downsampling 20 times.
## Let's make a simple calculation of the step lengths of real data:

tmp <- sampleMovement(real.data)

## and compute a good maximum allowed step length during optimization
## using the observed maximum divided by 20 (because each real step
## will comprise 20 simulated steps)

max.step.length <- max(tmp$stat[, "steplengths"]) / 20

## and finally build the species model with it.
## Note: "CRW.CRW.sl" is the short name for the model we want,
## as defined above

species.model <- speciesModel("CRW.CRW.sl", steplength = max.step.length)

## now run optimization

sol <- adjustModel(real.data, species.model, resol = 20, window = 10, nbins.hist = c(3, 3)
	, step.hist.log = TRUE)

## After finishing, we can extract the input parameters of the optimized solutions
## (100 by default) in the last generation (generation 1000 by default):

pars <- sol[[length(sol)]]$par

## now we can take the optimized solutions and reconstruct species based on them:

optimized.species <- apply(pars, 1, species.model)

## and make some simulations with those optimized species.
## Plot real trajectory

par(mfrow = c(2, 2), mar = c(0, 0, 1, 0))
plot(real.data, type = "l", asp = 1, axes = F, main = "Real")

## plot three simulated trajectories with optimized species

for(i in 1:3) {
	# remember we want to simulate at a 20 times higher frequency
	# so we do 500 (real data) x 20 steps
	sim <- simulate(optimized.species[[i]], 500 * 20)
	
	# now we downsample frequency to match real
	samp <- sampleMovement(sim, 20)
	
	# and plot the simulated trajectory before and after
	# downsampling 20 times.
	
	plot(sim[, 1:2], type = "l", asp = 1, axes = F, col = "gray"
		, main = "Simulated and downsampled")
	lines(samp$relocs, col = "black")
}

## Now plot the evolution of parameters along algorithm's generations.
## This is good to assess whether the final solutions converged
## but see ?generationPlot for details

generationPlot(sol, species.model)

}
}
\references{
  \itemize{
    \item{Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. 2002. A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation, 6(2), 182-197.}
  }
}

