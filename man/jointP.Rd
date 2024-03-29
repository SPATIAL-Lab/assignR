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
SpatRaster of probability density maps, e.g., as produced by \code{pdRaster}. All layers must have common spatial extent and projection.
}
}

\value{
SpatRaster.
}

\examples{
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

# joint probability for individuals of common origin
jointP(pp)
}
