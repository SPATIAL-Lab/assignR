\name{ham}

\alias{ham}

\docType{data}

\title{
data: adjacency matrix for keratin H isotope standards
}

\description{
This matrix documents keratin H isotope standard scales between which comparisons can be drawn, either through use of common reference materials or through comparative analyses.
}

\usage{data("ham")
}

\format{
matrix.
\describe{
  n x n symmetric, where n is the number of standard scales represented here and in \code{\link{hsds}}.
}
}

\references{
Magozzi, et al. (in review) A global compilation of known-origin keratin hydrogen and oxygen isotope data for wildlife and forensic applications.
}

\seealso{
\code{\link{hsds}}
}

\examples{
data("ham")
row.names(ham)
}

\keyword{datasets}
