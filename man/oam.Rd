\name{oam}

\alias{oam}

\docType{data}

\title{
data: adjacency matrix for keratin O isotope calibrations
}

\description{
This matrix documents keratin O isotope calibrations between which comparisons can be drawn, either through use of common reference materials or through comparative analyses.
}

\usage{data("oam")
}

\format{
matrix.
\describe{
  n x n symmetric, where n is the number of calibrations represented here and in \code{\link{orms}}.
}
}

\references{
Magozzi, et al. (in review) Calibration hierarchy transformation to improve the comparability of tissue organic hydrogen and oxygen isotope data.
}

\seealso{
\code{\link{orms}}
}

\examples{
data("oam")
row.names(oam)
}

\keyword{datasets}
