\name{knownOrig_sites}

\alias{knownOrig_sites}

\docType{data}

\title{
data: site data for tissues of known origin
}

\description{
This dataset consists of site locations and metadata for loactions were human hair, insect wings, and bird feathers of known origin have been collected.
}

\usage{data("knownOrig_sites")
}

\format{
SpatialPointsDataFrame with 5 fields. WGS84 unprojected geometry.
\describe{
  \item{[, 1]}{Site_ID: Unique ID}
  \item{[, 2]}{Site_name: Site name or descriptor}
  \item{[, 3]}{State: State or province of collection site, where recorded}
  \item{[, 4]}{Country: Country of collection site, where recorded}
  \item{[, 5]}{Site_comments: Site comments}
}
}

\references{
See \code{\link[assignR]{knownOrig_sources}} for citations for all datasets.
}

\seealso{
\code{\link[assignR]{knownOrig_samples}}, \code{\link[assignR]{knownOrig_sources}}, \code{\link[assignR]{subOrigData}}
}

\examples{
library(sp)
data("knownOrig_sites")
class(knownOrig_sites)
summary(knownOrig_sites)
plot(assignR:::wrld_simpl, border = "grey")
plot(knownOrig_sites, add = TRUE)
}

\keyword{datasets}
