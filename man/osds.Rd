\name{osds}

\alias{osds}

\docType{data}

\title{
Standards description for keratin O isotope standard scales
}

\description{
This data object contains information on keratin O isotope standard scales used across multiple laboratories since the year 2000.
}

\usage{data("osds")
}

\format{
data.frame with 17 fields.
\describe{
  \item{[, 1]}{Scale: Standard scale code}
  \item{[, 2]}{High_ID: Identifier for high-value standard}
  \item{[, 3]}{High_material: Description of high-value standard material}
  \item{[, 4]}{High_d18O: Mean oxygen isotope value of high-value standard}
  \item{[, 5]}{High_d18O_sd: Standard deviation of calibration data for high-value standard}
  \item{[, 6]}{High_d18O_n: Number of calibration data for high-value standard}
  \item{[, 7]}{High_d18O_se: Standard error of the calibrated mean for high-value standard}
  \item{[, 8]}{Low_ID: Identifier for low-value standard}
  \item{[, 9]}{Low_material: Description of low-value standard material}
  \item{[, 10]}{Low_d18O: Mean oxygen isotope value of low-value standard}
  \item{[, 11]}{Low_d18O_sd: Standard deviation of calibration data for low-value standard}
  \item{[, 12]}{Low_d18O_n: Number of calibration data for low-value standard}
  \item{[, 13]}{Low_d18O_se: Standard error of the calibrated mean for low-value standard}
  \item{[, 14]}{Citation: }
  \item{[, 15]}{Reference: }
  \item{[, 16]}{Ref_scale: Standard scale against which the values for this scale are calibrated}
  \item{[, 17]}{Treatment: Description of calibration proceedure}
}
}

\references{
Magozzi, et al. (in review) A global compilation of known-origin keratin hydrogen and oxygen isotope data for wildlife and forensic applications.
}

\seealso{
\code{\link{oam}}
}

\examples{
data("osds")
print(osds)
}

\keyword{datasets}
