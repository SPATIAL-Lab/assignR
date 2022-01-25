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
RasterStack or RasterBrick of probability density maps, e.g., as produced by \code{pdRaster}. All layers must have common spatial extent and projection.
}
}

\examples{
# load North America boundary and isoscape
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

# joint probability for individuals of common origin
jointP(pp)
}
