\name{jointP}

\alias{jointP}

\title{
Joint probability of origin
}

\description{
Joint probability for individuals of common origin (product of probabilities)
}

\usage{
jointP(pdR)
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
d = subOrigData(group = "Modern human", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap)

# joint probability for individuals of common origin
jointP(asn)
}
