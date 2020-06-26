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
data.frame with 18 fields.
\describe{
  \item{[, 1]}{Scale: Standard scale code}
  \item{[, 2]}{High_ID: Identifier for high-value standard}
  \item{[, 3]}{High_material: Description of high-value standard material}
  \item{[, 4]}{High: Mean oxygen isotope value of high-value standard}
  \item{[, 5]}{High_sd: Standard deviation of calibration data for high-value standard}
  \item{[, 6]}{High_n: Number of calibration data for high-value standard}
  \item{[, 7]}{High_se: Standard error of the calibrated mean for high-value standard}
  \item{[, 8]}{Low_ID: Identifier for low-value standard}
  \item{[, 9]}{Low_material: Description of low-value standard material}
  \item{[, 10]}{Low: Mean oxygen isotope value of low-value standard}
  \item{[, 11]}{Low_sd: Standard deviation of calibration data for low-value standard}
  \item{[, 12]}{Low_n: Number of calibration data for low-value standard}
  \item{[, 13]}{Low_se: Standard error of the calibrated mean for low-value standard}
  \item{[, 14]}{Ref_scale: Standard scale against which the values for this scale are calibrated}
  \item{[, 15]}{Citation_val: Source for the calibrated values used for this scale}
  \item{[, 16]}{Citation_cal: Source for the calibration methodology for this scale}
  \item{[, 17]}{Treatment: Description of calibration procedure}
  \item{[, 18]}{O_scale_comments: Scale comments}
}
}

\references{
Magozzi, S. et al. (in prep) Standard-chain transformation to improve the comparability of tissue hydrogen and oxygen isotope data.
}

\seealso{
\code{\link{oam}}
}

\examples{
data("osds")
print(osds)
}

\keyword{datasets}
