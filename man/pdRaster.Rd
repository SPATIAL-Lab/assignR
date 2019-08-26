\name{pdRaster}
\alias{pdRaster}
\title{
assign function
}
\description{
Assign unknown tissue isotope to calibrated precipitation isoscape.
}
\usage{
pdRaster(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, saveFiles = FALSE)
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
\item{saveFiles}{
logical. Save the results?
}
}
\value{
returns an object of class RasterStack or RasterLayer that includes a probability density raster for each individual in "unknown". If saveFiles == TRUE writes individual rasters in GTiff format to disk in the working directory for each individual in "unknown" and a pdf file with images for each probability density raster.
}

\seealso{
\code{\link[assignR]{calRaster}}
}
\examples{
data(naMap)
data(d2h_world)

#### Load hydrogen isotope for human hair in North America ####
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

#### Sample to assign
id <- "smile"
d2H <- -80
un <- data.frame(id,d2H)

asn <- pdRaster(r, unknown = un, mask = naMap, saveFiles = TRUE)
}
