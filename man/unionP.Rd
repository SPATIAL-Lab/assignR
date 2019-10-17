\name{unionP}
\alias{unionP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Union probability of origin
}
\description{
Probability that at least one individual came from the location (union of probabilities)
}
\usage{
unionP(pdR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdR}{
probability density map which must be raster or RasterLayer
}
}

\examples{
# load North America SpatialPolygonsDataFrame data in assignR package
data("naMap")
data("d2h_world")

#### Load hydrogen isotope for human hair in North America ####
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap)

# joint probability for individuals of common origin
unionP(asn)
}
