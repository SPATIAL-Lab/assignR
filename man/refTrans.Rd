\name{refTrans}

\alias{refTrans}

\title{
Transform calibration scale of data
}

\description{
This function conducts transformations to convert isotope measurements between calibration scales.
}

\usage{
refTrans(samples, marker = "d2H", std_scale = "VSMOW_H", niter = 1000)
}

\arguments{
  \item{samples}{
  data.frame. Must include a fields with data to be transformed, analtyical reproducibility of sample data (1 standard deviation), and original calibration scale for data. These fields must be named \emph{marker}, \emph{marker.sd}, and \emph{marker_std_scale}, respectively, where \emph{marker} is \dQuote{d2H} or \dQuote{d18O}. Values for the std_scale field should correspond to scale codes found in \code{\link[assignR]{hsds}} or \code{\link[assignR]{osds}}.
}
  \item{marker}{
character string. Column name for isotopic data to be extracted, either \dQuote{d2H} or \dQuote{d18O}.
}
  \item{std_scale}{
character string. Text identifier for standard calibration scale to which all isotope values will be transformed. See \code{\link[assignR]{hsds}} and \code{\link[assignR]{osds}}.}
  \item{niter}{
integer. Number of random samples used to propagate uncertainty in standard chain transformations.
}
}

\value{
list.
  \item{data}{
data.frame. Formatted identically to input object \emph{samples}, with values for the data and data uncertainty fields replaced with transformed values.}
  \item{chains}{
list. Each item is a character string vector containing the chain of standard scales used in the transformation for a set of selected samples. See \code{\link[assignR]{hsds}} and \code{\link[assignR]{osds}}.
}
}

\references{
Magozzi, S. et al. (in prep) Standard-chain transformation to improve the comparability of tissue hydrogen and oxygen isotope data.
}

\examples{
# Some fake sample data
s = data.frame("d2H" = seq(-100, -10, by=10), "d2H.sd" = rep(2), "d2H_std_scale" = rep("OldUT_H_1"))

# Transform to VSMOW-SLAP scale using default arguments
d1 = refTrans(s)

# Transformed values
d1$data$d2H

# error - target scale not valid for marker
\dontrun{d2 = refTrans(s, std_scale = "VSMOW_O")
}
}
