\name{knownOrig}
\alias{knownOrig}
\docType{data}
\title{
Isotope values for tissues of known origin
}
\description{
This dataset consists of hydrogen isotope values of human hair, insect chitin and bird feathers with location information.
}
\usage{data("knownOrig")}
\format{
  A SpatialPointsDataFrame
}
\seealso{
\code{\link[assignR]{subOrigData}}
}
\examples{
library(maptools)
data("knownOrig")
data("wrld_simpl")
plot(wrld_simpl)
points(knownOrig$Longitude,knownOrig$Latitude,col="red",cex=0.5)
}
\keyword{datasets}
