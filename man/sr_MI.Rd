\name{sr_MI}

\alias{sr_MI}

\docType{data}

\title{
data: low resolution locally weathered Sr isoscape for Michigan
}

\description{
Modeled 87Sr/86Sr value of the local rock weathering flux, obtained from waterisotopes.org and aggregated to 10 km resolution.
}

\usage{data("sr_MI")
}

\format{
RasterBrick with two layers. The first layer is the mean prediction and the second is 1 standard deviation (here estimated as 1 percent of the modeled mean)
}

\references{
Bataille, C. P. and Bowen, G. J. (2012) Mapping 87Sr/86Sr variations in bedrock and water for large scale provenance studies. \emph{Chemical Geology}, \bold{304--305}, 39--52.
}

\examples{
library(raster)

data("sr_MI")
plot(sr_MI)
}

