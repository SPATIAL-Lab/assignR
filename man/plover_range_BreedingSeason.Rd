\name{plover_range_BreedingSeason}

\alias{plover_range_BreedingSeason}

\docType{data}

\title{
Mountain plover breeding range}

\description{
Geographic range of North American mountain plover (\emph{Charadrius montanus}) during the breeding season.
}

\usage{data("plover_range_BreedingSeason")}

\format{SpatialPolygonsDataFrame}

\source{
\url{http://explorer.natureserve.org/servlet/NatureServe?searchName=Charadrius+montanus}
}

\references{
Ridgely, R. S. et. al (2003) Digital Distribution Maps of the Birds of the Western Hemisphere, version 1.0. NatureServe, Arlington, Virginia, USA.
}

\examples{
data("plover_range_BreedingSeason")
plot(plover_range_BreedingSeason)
}

\keyword{datasets}
