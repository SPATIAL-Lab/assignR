\name{ham}

\alias{ham}

\docType{data}

\title{
data: adjacency matrix for keratin H isotope calibrations
}

\description{
This matrix documents keratin H isotope calibrations between which comparisons can be drawn, either through use of common reference materials or through comparative analyses
}

\usage{data("ham")
}

\format{
matrix.
\describe{
  n x n symmetric, where n is the number of calibrations represented here and in \code{\link{hrms}}.
}
}

\references{
Magozzi, et al. (in review) Calibration hierarchy transformation to improve the comparability of organic hydrogen and oxygen isotope data.
}

\seealso{
\code{\link{hrms}}
}

\examples{
data("ham")
row.names(ham)
}

\keyword{datasets}
