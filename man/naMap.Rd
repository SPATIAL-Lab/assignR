\name{naMap}
\alias{naMap}
\docType{data}
\title{
data: naMap
}
\description{
North America polygon
}
\usage{data("naMap")}
\format{
  The format is SpatialPolygonsDataFrame, with extent of -164.7, -51.7, 21, 73  (xmin, xmax, ymin, ymax). coord. ref. is "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

}
\examples{
data("naMap")
plot(naMap)
}
