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
RasterStack or RasterBrick of probability density maps, e.g., as produced by \code{pdRaster}.
}
  \item{sites}{
SpatialPoints* object containing the collection locations for the samples represented in pdR.
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


}

