\name{unionP}

\alias{unionP}

\title{
Union probability of origin
}

\description{
Probabilities that at least one individual came from each location in the assignment area (union of probabilities)
}

\usage{
unionP(pdR)
}

\arguments{
  \item{pdR}{
SpatRaster of probability density maps, e.g., as produced by \code{pdRaster}.
}
}

\value{
SpatRaster.
}

\examples{
library(terra)

# load North America boundary and global isoscape
data("naMap")
d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, genplot = FALSE)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(d, d2h_lrNA, naMap, genplot = FALSE)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id, d2H)

# posterior probabilities
pp = pdRaster(r, un, mask = naMap, genplot = FALSE)

# probability that one or more individuals are from a given location
unionP(pp)
}
