\name{knownOrig}

\alias{knownOrig}

\docType{data}

\title{
data: database of H and O isotope data for tissues of known origin
}

\description{
This dataset consists of hydrogen and oxygen isotope values and metadata for human hair, insect wings, and bird feathers of known geographic origin.
}

\usage{data("knownOrig")
}

\format{
list.
\describe{
\item{sites}{SpatialPointsDataFrame with 5 fields. WGS84 unprojected geometry.
  \describe{
    \item{[, 1]}{Site_ID: Unique ID}
    \item{[, 2]}{Site_name: Site name or descriptor}
    \item{[, 3]}{State: State or province of collection site, where recorded}
    \item{[, 4]}{Country: Country of collection site, where recorded}
    \item{[, 5]}{Site_comments: Site comments}
  }
}
\item{samples}{data.frame with 15 fields.
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
\item{sources}{data.frame with 17 fields.
  \describe{
    \item{[, 1]}{Dataset_ID: Unique ID}
    \item{[, 2]}{Dataset_name: Short name or descriptor}
    \item{[, 3]}{Citation: Bibliographic citation for study}
    \item{[, 4]}{Sampling_method: How material was subsampled for analysis, if reported}
    \item{[, 5]}{Sample_powdered: Was sample powdered prior to analysis (Y/N/NA)?}
    \item{[, 6]}{Lipid_extraction: Were lipids chemically extracted prior to analysis (Y/N/NA)?}
    \item{[, 7]}{Lipid_extraction_method: Solvent used to extract lipids}
    \item{[, 8]}{Exchange: Was a correction for exchangeable H made (Y/N/NA)?}
    \item{[, 9]}{Exchange_method: Method used to correct for exchangeable H}
    \item{[, 10]}{Exchange_T: Was H exchange carried out at ambient or high temperature (Ambient/High/NA)?}
    \item{[, 11]}{H_cal: Reference scale used to calibrate H isotope data, see \code{\link[assignR]{stds}} object hstds}
    \item{[, 12]}{O_cal: Reference scale used to calibrate O isotope data, see \code{\link[assignR]{stds}} object ostds}
    \item{[, 13]}{Std_powdered: Were calibration standards powdered (Y/N/NA)?}
    \item{[, 14]}{Drying: Did the study document how samples were fully dried and transferred dry to instrument (Y/N/NA)?}
    \item{[, 15]}{Analysis_method: Instrument configuration used for analysis}
    \item{[, 16]}{Analysis_type: What elements were analyzed for stable isotope ratios (H/O/H_O)?}
    \item{[, 17]}{Source_comments: Data source comments}
  }
}
}
}

\seealso{
\code{\link[assignR]{subOrigData}}
}

\examples{
library(sp)
data("knownOrig")
class(knownOrig$sites)
class(knownOrig$samples); class(knownOrig$sources)
summary(knownOrig$samples)
print(knownOrig$sources[, 1:2])
plot(assignR:::wrld_simpl, border = "grey")
plot(knownOrig$sites, add = TRUE)
}

\keyword{datasets}