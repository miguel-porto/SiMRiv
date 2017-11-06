\name{binCounts}
\alias{binCounts}
\title{Count values in given bins}
\description{
  A convenience function to count the values (optionally log-transformed) falling within given equal-range bins.
}
\usage{
binCounts(data, range, nbins, log = FALSE)
}
\arguments{
  \item{data}{a numeric vector with the data.}
  \item{range}{the closed data range in which to divide bins, as a two-element numeric vector. Values falling outside \code{range} are discarded.}
  \item{nbins}{the number of bins to split the data. Result is guaranteed to have this number of bins, irrespective of data.}
  \item{log}{whether or not to log-transform data and ranges before binning.}
}
\details{
  This is just a convenience function to simplify the task of couting data in equal bins.
  The same result could be achieved with \code{hist} with the right arguments.
  It is mostly used internally by \code{\link{adjustModel}} during the approximation on input parameters,
  for the computation of objectives.
}
\value{
	A named vector with the counts of values in each bin.
}
\seealso{
  \code{\link{adjustModel}}.
}

