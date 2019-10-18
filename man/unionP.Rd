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
RasterStack or RasterBrick of probability density maps, e.g., as produced by \code{calRaster}. All layers must have common spatial extent and projection.
}
}

\examples{
# load North America boundary and global isoscape
data("naMap")
data("d2h_world")

# load hydrogen isotope data for human hair in North America
d = subOrigData(taxon = "Homo sapiens", 
  reference = "Ehleringer et al. 2008", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap)

# probability that one or more individuals are from a given location
unionP(asn)
}
