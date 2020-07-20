\name{oam}

\alias{oam}

\docType{data}

\title{
data: adjacency matrix for keratin O isotope standards
}

\description{
This matrix documents keratin O isotope standard scales between which comparisons can be drawn, either through use of common reference materials or through comparative analyses.
}

\usage{data("oam")
}

\format{
matrix.
\describe{
  n x n symmetric, where n is the number of standard scales represented here and in \code{\link{osds}}.
}
}

\references{
Magozzi, et al. (in review) A global compilation of known-origin keratin hydrogen and oxygen isotope data for wildlife and forensic applications.
}

\seealso{
\code{\link{osds}}
}

\examples{
data("oam")
row.names(oam)
}

\keyword{datasets}
