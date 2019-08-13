\name{d18o_world}
\alias{d18o_world}
\docType{data}
\title{
data: d18o_world
}
\description{
This file is the global growing season d18o precipitation isoscape available on waterisotopes.org and reported in Bowen et al., 2005.
}
\usage{data("d18o_world")}
\format{
  raster with two layers. The first layer is the mean and the second layer is the 2 standard deviation
}

\references{
Bowen, G.J., Wassenaar, L.I. and Hobson, K.A., 2005. Global application of stable hydrogen and oxygen isotopes to wildlife forensics. Oecologia, 143(3), pp.337-348.
}
\examples{
data(d18o_world)
plot(d18o_world)
}
