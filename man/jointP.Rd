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
probability density map which must be raster or RasterLayer
}
}

\examples{
# load North America SpatialPolygonsDataFrame data in assignR package
data("naMap")
data(d2h_world)
d = subOrigData(taxon = c("Homo sapiens"), mask = naMap)

# exclude some outlier. This step is optional, which depends on your data quality.
dd = d[d@coords[,1] < -80,]

# Rescale from environmental isoscape to tissue isoscape
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)

# four unknown-origin example
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap, saveFile = FALSE)

# joint probability for individuals of common origin
jointP(asn)
}
