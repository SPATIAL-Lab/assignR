\name{pdRaster}

\alias{pdRaster}

\title{
Probability of origin surfaces
}

\description{
Calculate posterior probabilities of origin for a sample based on its isotope ratio.
}

\usage{
pdRaster(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, outDir = NULL)
}

\arguments{
\item{r}{
SpatRaster with two layers, \code{rescale} object (see \code{\link{calRaster}}), or \code{\link{isoStack}} object. For user-generated raster objects, the first layer must be the substrate-specific isoscape (mean prediction) and the second the isoscape prediction uncertainty (1 standard deviation).
}
\item{unknown}{
data.frame, \code{\link{refTrans}} object, or list of two or more \code{\link{refTrans}} objects. For user-created data.frame, first column should contain unique sample IDs, and subsequent columns should contain sample isotope values for one or more isotopes to be used in assignment.
}
\item{prior}{
SpatRaster. Optional raster layer with prior probabilities, which has the same projection, resolution and extent as \code{r}.
}
\item{mask}{
SpatialPolygonsDataFrame. This mask will constrain the assignment area. If this is not provided, a default of mask of the extent of \code{r} is used.
}
\item{genplot}{
logical. Plot results in R.
}
\item{outDir}{
character string. Directory path to which output will be saved. If NULL no files are written.
}
}

\details{
If more than one isotope marker is to be used for multivariate assignment, \code{r} must be an \code{\link{isoStack}} object and the number of isoscapes in that object must be equal to the number of isotope-value columns or \code{\link{refTrans}} objects included in \code{unknown}. Isoscapes and unknown sample values will be matched based on order, so it is critical that the values appear in the same order in these two input objects.
}

\value{
SpatRaster including a probability density surface for each individual in \code{unknown}. If \code{outDir} is not NULL, writes individual rasters in GeoTIFF format and a single PDF file with images for each probability density raster to the designated directory.
}

\seealso{
\code{\link{calRaster}}
\code{\link{refTrans}}
\code{\link{isoStack}}
}

\examples{
library(terra)

# load North America boundary and global isoscape
data("naMap")
d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(d, d2h_lrNA, naMap)

# sample to assign
id = "smile"
d2H = -80
un = data.frame(id, d2H)

# posterior probability surface
asn = pdRaster(r, un, mask = naMap)
}
