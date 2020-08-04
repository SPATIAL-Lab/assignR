\name{knownOrig_sources}

\alias{knownOrig_sources}

\docType{data}

\title{
data: data source information  for tissues of known origin
}

\description{
This dataset consists of citation, methodological, and other metadata for studies from which data on known origin human hair, insect wings, and bird feathers have been compiled.
}

\usage{data("knownOrig_sources")
}

\format{
data.frame with 17 fields.
\describe{
  \item{[, 1]}{Dataset_ID: Unique ID}
  \item{[, 2]}{Dataset_name: Dataset name or descriptor}
  \item{[, 3]}{Citation: Bibliographic citation for study}
  \item{[, 4]}{Sampling_method: How material was subsampled for analysis, if reported}
  \item{[, 5]}{Sample_powdered: Was sample powdered prior to analysis (Y/N/NA)?}
  \item{[, 6]}{Lipid_extraction: Were lipids chemically extracted prior to analysis (Y/N/NA)?}
  \item{[, 7]}{Lipid_extraction_method: Solvent used to extract lipids}
  \item{[, 8]}{Exchange: Was a correction for exchangeable H made (Y/N/NA)?}
  \item{[, 9]}{Exchange_method: Method used to correct for exchangeable H}
  \item{[, 10]}{Exchange_T: Was H exchange carried out at ambient or high temperature (Ambient/High/NA)?}
  \item{[, 11]}{H_cal: Reference scale used to calibrate H isotope data, see \code{\link[assignR]{hrms}}}
  \item{[, 12]}{O_cal: Reference scale used to calibrate O isotope data, see \code{\link[assignR]{orms}}}
  \item{[, 13]}{Std_powdered: Were calibration standards powdered (Y/N/NA)?}
  \item{[, 14]}{Drying: Did the study document that samples were fully dried and transferred dry to instrument (Y/N/NA)?}
  \item{[, 15]}{Analysis_method: Instrument configuration used for analysis}
  \item{[, 16]}{Analysis_type: What elements were analyzed for stable isotope ratios (H/O/H_O)?}
  \item{[, 17]}{Source_comments: Data source comments}
}
}

\seealso{
\code{\link[assignR]{knownOrig_samples}}, 
\code{\link[assignR]{knownOrig_sites}}, 
\code{\link[assignR]{subOrigData}}
}

\examples{
data("knownOrig_sources")
class(knownOrig_sources)
summary(knownOrig_sources)
}

\keyword{datasets}
