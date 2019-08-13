\name{jointP}
\alias{jointP}
\title{
jointP
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
# Rescale from environmental isoscape to tissue isoscape (known must be SpatialPointsDataFrame with coordinate reference system )
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)
# four unknown-origin example
id = c("A", "B", "C", "D")
d2H = c(-110, -180, -130, -150)
un = data.frame(id,d2H)
# assignment
asn = pdRaster(r,unknown=un,mask=naMap)
# joint probability for individuals of common origin
jointP(asn)
}
