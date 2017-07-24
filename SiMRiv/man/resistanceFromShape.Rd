\name{resistanceFromShape}
\alias{resistanceFromShape}
\title{Build resistance raster by combining shapefiles}
\description{
  Creates a resistance raster to be used in simulations, by rasterizing and combining different shapefiles.
  It is basically a helper function that uses the functions from package \code{\link{raster-package}} to create and manipulate
  such raster.
}
\usage{
  resistanceFromShape(shp, baseRaster, res, binary = is.na(field)
	, field = NA, background = 1, buffer = NA, margin = 0
	, mapvalues = NA, extend = TRUE, ...)
}
\arguments{
  \item{shp}{either a character string specifying the shapefile filename or a shapefile object itself.}
  \item{baseRaster}{if provided, a raster onto which to stack the given rasterized shapefile. If not provided, a new raster will be created.}
  \item{res}{the desired pixel resolution of the raster to be created, when \code{baseRaster} is not provided.}
  \item{binary}{if TRUE, the shapefile will be rasterized so that any feature is assigned a value of 1, and the background 0.}
  \item{field}{either a number in the range [0-1], in which case it will be assigned to all pixels covered by features of the shapefile; or a character string naming the shapefile field from which to extract such number.}
  \item{background}{the value in the range [0-1] to assign to all pixels that are not covered by any shapefile feature.}
  \item{buffer}{the size of a buffer to build around every shapefile feature before rasterizing.}
  \item{margin}{the margin to leave around the shapefile's extent when rasterizing (i.e. how much to increase shapefile's extent).}
  \item{mapvalues}{a named vector specifying the resistance value mapping from the classes of \code{field}.}
  \item{extend}{set to TRUE to extend \code{baseRaster} if the shapefile has a larger extent. If FALSE, the shapefile will be clipped to \code{baseRaster}'s extent.}
  \item{...}{other arguments to pass to \code{\link{rasterize}}.}
}
\details{
  This function rasterizes the given shapefile using provided options and optionally stacks it onto the provided \code{baseRaster}.
}
\note{
  This function
}
\value{
  The resistance raster, an object of class \code{\link{RasterLayer-class}}.
}
\seealso{
  \code{\link{simulate}}.
}
\examples{

}
\references{
  \itemize{
  }
}
