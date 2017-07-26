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
  \item{binary}{if TRUE, the shapefile will be rasterized so that any feature is assigned a value of 0, and the background 1.}
  \item{field}{either a number in the range [0-1], in which case it will be assigned to all pixels covered by features of the shapefile; or the name of the numeric shapefile field from which to extract such number; or the name of the factor shapefile field containing classes to map to resistance values (see \code{mapvalues}).}
  \item{background}{the value in the range [0-1] to assign to all pixels that are not covered by any shapefile feature.}
  \item{buffer}{the size of a buffer to build around every shapefile feature before rasterizing.}
  \item{margin}{the margin to leave around the shapefile's extent when rasterizing (i.e. how much to increase shapefile's extent).}
  \item{mapvalues}{a named vector specifying the resistance value mapping from the classes of \code{field}.}
  \item{extend}{set to TRUE to extend \code{baseRaster} if the shapefile has a larger extent. If FALSE, the shapefile will be clipped to \code{baseRaster}'s extent.}
  \item{...}{other arguments to pass to \code{\link{rasterize}}.}
}
\details{
  This function rasterizes the given shapefile using provided options and optionally stacks it onto the provided \code{baseRaster}.
  The produced raster does not contain NAs and all values are in the range [0, 1]. All the areas of the raster for which data is not provided are assigned the value of \code{background}.
  
  When combining a shapefile to an existing \code{baseRaster}, only the areas covered by features are updated in the base raster;
  all the remaining areas are left with the original values of \code{baseRaster}. If the shapefile to combine has a larger extent than \code{baseRaster},
  those extra pixels will be assigned the \code{background} value defined for the shapefile (not use the original background of \code{baseRaster}).
}
\value{
  The resistance raster, an object of class \code{\link{RasterLayer-class}}.
}
\seealso{
  \code{\link{simulate}}.
}
\examples{
## Example taken from the vignette; see the vignette
## for more details and examples.
## In this example we read a land cover shapefile and
## assign resistance values based on each polygon's
## land cover class (provided in the field 'coverclass')

landcover <- resistanceFromShape(
  system.file("doc/landcover.shp", package="SiMRiv")
  , res = 150, field = "coverclass", mapvalues = c(
    "forest" = 0.5, "urban" = 1, "dam" = 0
    , "shrubland" = 0.75), background = 0.95)
    
## We then combine it with a river network from another
## shapefile, assigning a value of 0 to all rivers

river.landcover <- resistanceFromShape(
  system.file("doc/river-sample.shp", package="SiMRiv")
  , baseRaster = landcover, buffer = 100, field = 0
  , background = 0.95, margin = 1000)

plot(river.landcover, axes = FALSE, mar = c(0, 0, 0, 2))

}

