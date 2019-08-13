\name{d2h_world}
\alias{d2h_world}
\docType{data}
\title{
data: d2h_world
}
\description{
precipitation δ2H isoscape of global growing season from waterisotopes.org
}
\usage{data("d2h_world")}
\format{
 raster with two layers. The first layer is mean and the second one is 2 standard deviation
}
\references{
Bowen, G.J., Wassenaar, L.I. and Hobson, K.A., 2005. Global application of stable hydrogen and oxygen isotopes to wildlife forensics. Oecologia, 143(3), pp.337-348.
}
\examples{
data(d2h_world)
plot(d2h_world)
}

