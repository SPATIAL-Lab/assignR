\name{wDist}

\alias{wDist}

\title{
Probability weighted distances and bearings
}

\description{
Calculate the distance and bearing of migration for one or more samples, weighted by probabilities from a \code{pdRaster} analysis. 
}

\usage{
wDist(pdR, sites)
}

\arguments{
  \item{pdR}{
RasterStack or RasterBrick of n probability density maps, e.g., as produced by \code{pdRaster}.
}
  \item{sites}{
SpatialPoints* object containing the collection locations for the n samples represented in pdR.
}
}

\details{
\code{pdR} and \code{sites} must be of equal length and corresponding order. Names in the output object are taken from the names of the layers in \code{pdR}.

Distances and bearings are calculated on the WGS84 geoid using functions from the \pkg{geosphere} package.

Bearing values correspond to the intial bearing from source to collection location, and are reported on a scale of -180 to +180 degrees. The statistical metrics are rectified so that values for distributions spanning due south are reported correctly. Both weighted bearing and distance distributions are often multimodal, and it is recommended to review the distribution densities to assess the representativeness of the statitics (e.g., using \code{\link{plot.wDist}}). 
}

\value{
Returns an object of class \dQuote{wDist}, a list of length n.
  \item{stats}{
named number. Statistics for the distance (dist, meters) and bearing (bear, degrees) between source and collection locations, including the weighted mean (wMean) and quantile (wXX) values.}
  \item{d.density}{
density. Weighted kernal density for the distance between source and collection locations (meters). See \code{\link{stats::density()}}. 
  }
  \item{s.density}{
density. Weighted kernal density for the bearing between source and collection locations (degrees). See \code{\link{stats::density()}}.
  }
}

\examples{
library(sp)

# load North America boundary and global isoscape
data("naMap")
data("d2h_lrNA")

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, niter = 100)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# posterior probabilities
pp = pdRaster(r, unknown = un, mask = naMap)

# random collection locations
sites = d$data[sample(seq(length(d$data)), 5),]

# generate a wDist object
wd = wDist(pp, sites)

# combine stats and print
c(wd)
}

