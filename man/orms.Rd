\name{orms}

\alias{orms}

\docType{data}

\title{
data: description for keratin H isotope reference material calibrations
}

\description{
This data object contains information on keratin O isotope reference materials and RM calibrations used across multiple laboratories since the year 2000.
}

\usage{data("orms")
}

\format{
data.frame with 18 fields.
\describe{
  \item{[, 1]}{Calibration: Calibration code}
  \item{[, 2]}{High_ID: Identifier for high-value RM}
  \item{[, 3]}{High_material: Description of high-value RM material}
  \item{[, 4]}{High: Mean oxygen isotope value of high-value RM}
  \item{[, 5]}{High_sd: Standard deviation of calibration data for high-value RM}
  \item{[, 6]}{High_n: Number of calibration data for high-value RM}
  \item{[, 7]}{High_se: Standard error of the calibrated mean for high-value RM}
  \item{[, 8]}{Low_ID: Identifier for low-value RM}
  \item{[, 9]}{Low_material: Description of low-value RM material}
  \item{[, 10]}{Low: Mean oxygen isotope value of low-value RM}
  \item{[, 11]}{Low_sd: Standard deviation of calibration data for low-value RM}
  \item{[, 12]}{Low_n: Number of calibration data for low-value RM}
  \item{[, 13]}{Low_se: Standard error of the calibrated mean for low-value RM}
  \item{[, 14]}{Ref_scale: Calibration scale against which the values for this calibration are anchored}
  \item{[, 15]}{Citation_val: Source for the calibrated values}
  \item{[, 16]}{Citation_cal: Source for the methodology used for this calibration}
  \item{[, 17]}{Treatment: Description of calibration procedure}
  \item{[, 18]}{O_calibration_comments: Comments}
}
}

\references{
Magozzi, S. et al. (in prep) Calibration hierarchy transformation to improve the comparability of organic hydrogen and oxygen isotope data.
}

\seealso{
\code{\link{oam}}
}

\examples{
data("orms")
print(orms)
}

\keyword{datasets}
