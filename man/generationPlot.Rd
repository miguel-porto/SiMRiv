\name{generationPlot}
\alias{generationPlot}
\title{Plots input parameter optimization results}
\description{
  Plots the evolution of the optimized solutions (sets of input parameters) along the \code{\link{adjustModel}} algorithm's generations.
}
\usage{
  generationPlot(solutions, species.model
    , plot.quantiles = c(0.10, 0.5, 0.90), only.pareto = FALSE
    , show.legend = TRUE, lwd = 1.5, mar = c(2.3, 2.3, 0.2, 2.3)
    , mgp = c(1.2, 0.2, 0), tcl = -0.25, ...)
}
\arguments{
  \item{solutions}{the result of \code{\link{adjustModel}}.}
  \item{species.model}{the species model that was adjusted and that was passed to \code{\link{adjustModel}}.}
  \item{plot.quantiles}{the three quantiles to plot. The middle is plotted as a solid line, the extremes define the shaded area.}
  \item{only.pareto}{whether to display the quantiles only of the Pareto front (\code{TRUE}) or of all solutions (\code{FALSE}).}
  \item{show.legend}{whether to show a legend indicating the correspondence of colors - input parameters.}
  \item{lwd}{line width to be used to draw the middle quantile (usually the median).}
  \item{mar, mgp, tcl, ...}{arguments to pass to \code{\link{par}}.}
}
\details{
  The generation plot depicts, for each generation of the optimization algorithm, the given quantiles of each input parameter being optimized
  in the population of solutions. This is not the ideal plot because input parameters are plotted independently, while they are only supposed to
  make sense in the context of a given solution (i.e. it is the combination of the input parameters that is being optimized, not the parameters in isolation).
  
  Nonetheless, this plot is still a good way to assess whether the final solutions converged to stable values, which can indicate that the algorithm succeeded in replicating the real trajectory provided.
  However, note that depending on how the problem is formulated, the solutions might not be expected to converge to a single solution "type":
  for example, there can be two types of solutions, corresponding to two different ways of achieving similar results.
}
\value{
  Returns, invisibly, a 3-D matrix with the computed quantiles for all input parameters along generations.
}
\seealso{
  \code{\link{adjustModel}}.
}
\examples{
## see ?adjustModel for a complete example
}

