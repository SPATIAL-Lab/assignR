\name{knownOrig_samples}

\alias{knownOrig_samples}

\docType{data}

\title{
Sample data for tissues of known origin
}

\description{
This dataset consists of hydrogen and oxygen isotope values of human hair, insect wings, and bird feathers plus sample metadata.
}

\usage{data("knownOrig_samples")
}

\format{
data.frame with 15 fields.
\describe{
  \item{[, 1]}{Sample_ID: Unique ID}
  \item{[, 2]}{Sample_ID_orig: ID used in original data report}
  \item{[, 3]}{Site_ID: ID for sample collection site}
  \item{[, 4]}{Dataset_ID: ID for dataset from which sample is derived}
  \item{[, 5]}{Taxon: Genus and species name}
  \item{[, 6]}{Group: Biological group (informal)}
  \item{[, 7]}{Source_quality: Code indicating level of certainty in geographic origin}
  \item{[, 8]}{Age_class: Code for age of individual}
  \item{[, 9]}{Material_type: Tissue sampled, e.g., \dQuote{Hair}}
  \item{[, 10]}{Matrix: Compound measured, e.g., \dQuote{Keratin}}
  \item{[, 11]}{d2H: Hydrogen isotope value (permil)}
  \item{[, 12]}{d2H.sd: Reported analytical uncertainty for hydrogen isotope value (permil)}
  \item{[, 13]}{d18O: Oxygen isotope value (permil)}
  \item{[, 14]}{d18O.sd: Reported analytical uncertainty for oxygen isotope value (permil)}
  \item{[, 15]}{Sample_comments: Sample comments}
}
}

\references{
See \code{\link[assignR]{knownOrig_sources}} for citations for all datasets.
}

\seealso{
\code{\link[assignR]{knownOrig_sites}}
\code{\link[assignR]{knownOrig_sources}}
\code{\link[assignR]{subOrigData}}
}

\examples{
data("knownOrig_samples")
class(knownOrig_samples)
summary(knownOrig_samples)
}

\keyword{datasets}
