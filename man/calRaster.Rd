\name{calRaster}

\alias{calRaster}

\title{
Rescale isoscape using linear regression
}

\description{
This function uses known-origin tissue data to rescale a map of environmental isotope values to a map of tissue value (and associated uncertainty) using a linear regression model.
}

\usage{
calRaster(known, isoscape, mask = NULL, interpMethod = 2, NA.value = NA, 
  ignore.NA = TRUE, genplot = TRUE, outDir =  NULL, verboseLM = TRUE)
}

\arguments{
  \item{known}{
subOrigData or SpatialPointsDataFrame. Known-origin tissue isotope data from the subOrigData function or provided by user. User-provided data must be formatted as a subOrigData object (see \code{\link[assignR]{subOrigData}} or a SpatialPointsDataFrame in which the first data field contains the measured tissue isotope value and the second the 1 standard deviation uncertainty on that value.
}
  \item{isoscape}{
SpatRaster. Isoscape raster with two layers. The first one is the mean and the second is one standard deviation.
}
  \item{mask}{
SpatialPolygons or SpatialPolygonsDataFrame. Constrains the area of the output rasters. If this is not provided, the entire area of \code{isoscape} is returned.
}
  \item{interpMethod}{
numeric. 1 or 2. Designate one of two methods for extracting values from \code{isoscape}. If 1, values for the cell in which a point falls are returned. If 2 the returned values are interpolated (bilinear interpolation) from the values of the four adjacent grid cells.
}
\item{NA.value}{
NA or numeric. Value representing the absence of data in \code{isoscape}. Commonly used values include NA, -9999, and 9999.
}
\item{ignore.NA}{
logical. If NA values are extracted from \code{isoscape} at the location of samples in known, ignore these values and proceed with \code{calRaster}.
}
\item{genplot}{
logical. Plot the results.
}
\item{outDir}{
character string. Directory path to which output will be saved. If NULL no files are written.
}
\item{verboseLM}{
logical. Print out the linear regression results.
}
}

\value{
Returns an object of class \dQuote{rescale}.
\item{isoscape.rescale}{SpatRaster. \code{isoscape} rescaled using a linear regression model between values of \code{known} and \code{isoscape}. The first layer is the mean prediction and the second the standard deviation.}
\item{lm.data}{data.frame. Known origin data and extracted \code{isoscape} values used for linear regression modeling.}
\item{lm.model}{list. Linear regression model.}
}

\seealso{
\code{\link[assignR]{pdRaster}}
}

\examples{
# load North America boundary and global isoscape
library(terra)
data("naMap")
d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, niter = 100)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)
}
