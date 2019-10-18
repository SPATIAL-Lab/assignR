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
  ignore.NA = TRUE, genplot = TRUE, savePDF =  FALSE, verboseLM = TRUE)
}

\arguments{
  \item{known}{
SpatialPointsDataFrame. Known-origin data that should contain only one feature: tissue isotope value. Known-origin data can be queried using \code{knownOrig}.
}
  \item{isoscape}{
raster. Isoscape raster with two layers. The first one is the mean and the second is one standard deviation.
}
  \item{mask}{
SpatialPolygonsDataFrame. Constrains the area of the output rasters. If this is not provided, the entire area of \code{isoscape} is returned.
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
\item{savePDF}{
logical. Save the results as pdf.
}
\item{verboseLM}{
logical. Print out the linear regression results.
}
}

\value{
\item{isoscape.rescale}{RasterBrick. \code{isoscape} rescaled using a linear regression model between values of \code{known} and \code{isoscape}. The first layer is the mean prediction and the second the standard deviation.}
\item{lm.data}{data.frame. Known origin data and extracted \code{isoscape} values used for linear regression modeling.}
\item{lm.model}{list. Linear regression model.}
}

\seealso{
\code{\link[assignR]{pdRaster}}
}

\examples{
# load North America boundary and global isoscape
data("naMap")
data("d2h_world")

# load hydrogen isotope data for human hair in North America
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)
}
