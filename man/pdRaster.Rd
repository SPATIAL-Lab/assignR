\name{pdRaster}
\alias{pdRaster}
\title{
assign function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape.
}
\usage{
pdRaster(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, saveFile = TRUE)
}
\arguments{
\item{r}{
output of \code{\link[assignR]{calRaster}}, or RasterStack or RasterBrick with two layers. The tissue-calibrated isoscape raster (mean and sd) obtained from rescale function, or user provided
}
\item{unknown}{
data.frame. First column should be unique sample IDs, and second column should be tissue isotope values for the individuals for which the assignments will be made
}
\item{prior}{
raster. This is an optional raster layer with prior probabilities, which should has the same projection, resolution and extent as the isoscape raster
}
\item{mask}{
SpatialPolygonsDataFrame. This mask will constrain the assignment area. If this is not provided, a default of mask of rescaled.mean is used.
}
\item{genplot}{
logical. Generate summary plots in R?
}
\item{saveFile}{
logical. Save the results?
}
}
\value{
returns an object of class RasterStack or RasterLayer that includes a probability density raster for each individual in "unknown". pdRaster writes individual rasters in GTiff format to disk in the working directory for each individual in "unknown".
}

\seealso{
\code{\link[assignR]{calRaster}}
}
\examples{
data(naMap)
data(d2h_world)
d <- subOrigData(taxon = c("Homo sapiens"), mask = naMap)
r <- calRaster(known = d, isoscape = d2h_world, mask = naMap)
id <- "smile"
d2H <- -80
un <- data.frame(id,d2H)
asn <- pdRaster(r, unknown = un, mask=naMap)
}
