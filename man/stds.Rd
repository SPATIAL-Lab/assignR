\name{stds}

\alias{stds}

\docType{data}

\title{
data: documentation of keratin H and O isotope standard calibrations
}

\description{
This data object contains information on keratin H and O isotope standard materials and calibrations used across multiple laboratories since the year 2000.
}

\usage{
data("stds")
}

\format{
list.
\describe{
\item{hstds}{
  data.frame with 18 fields.
  \describe{
    \item{[, 1]}{Calibration: Calibration code}
    \item{[, 2]}{High_ID: Identifier for high-value standard}
    \item{[, 3]}{High_material: Description of high-value standard material}
    \item{[, 4]}{High: Mean hydrogen isotope value of high-value standard}
    \item{[, 5]}{High_sd: Standard deviation of calibration data for high-value standard}
    \item{[, 6]}{High_n: Number of calibration data for high-value standard}
    \item{[, 7]}{High_se: Standard error of the calibrated mean for high-value standard}
    \item{[, 8]}{Low_ID: Identifier for low-value standard}
    \item{[, 9]}{Low_material: Description of low-value standard material}
    \item{[, 10]}{Low: Mean hydrogen isotope value of low-value standard}
    \item{[, 11]}{Low_sd: Standard deviation of calibration data for low-value standard}
    \item{[, 12]}{Low_n: Number of calibration data for low-value standard}
    \item{[, 13]}{Low_se: Standard error of the calibrated mean for low-value standard}
    \item{[, 14]}{Ref_scale: Calibration scale against which the values for this calibration are anchored}
    \item{[, 15]}{Citation_val: Source for the calibrated values}
    \item{[, 16]}{Citation_cal: Source for the methodology used for this calibration}
    \item{[, 17]}{Treatment: Description of calibration procedure}
    \item{[, 18]}{H_calibration_comments: Comments}
  }
}
\item{ostds}{
  data.frame with 18 fields.
  \describe{
    \item{[, 1]}{Calibration: Calibration code}
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
    \item{[, 14]}{Ref_scale: Calibration scale against which the values for this calibration are anchored}
    \item{[, 15]}{Citation_val: Source for the calibrated values}
    \item{[, 16]}{Citation_cal: Source for the methodology used for this calibration}
    \item{[, 17]}{Treatment: Description of calibration procedure}
    \item{[, 18]}{O_calibration_comments: Comments}
  }
}
\item{ham}{
  matrix. n x n symmetric, where n is the number of calibrations represented here and in \code{stds$hstds}.
}
\item{oam}{
  matrix. n x n symmetric, where n is the number of calibrations represented here and in \code{stds$ostds}.
}
}
}

\references{
Magozzi, S. et al. (revised) Calibration chain transformation to improve the comparability of organic hydrogen and oxygen isotope data.
}

\examples{
library(raster)
data("stds")

print(stds$hstds[, 1:5])
print(stds$ostds[, 1:5])
image(stds$ham)
image(stds$oam)
}

\keyword{datasets}