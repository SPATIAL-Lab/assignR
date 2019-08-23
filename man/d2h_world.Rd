\name{d2h_world}
\alias{d2h_world}
\docType{data}
\title{
data: d2h_world
}
\description{
Global growing season precipitation δ2H isoscape from waterisotopes.org
}
\usage{data("d2h_world")}
\format{
 raster with two layers. The first layer is the mean prediction and the second one is 1 standard deviation
}
\references{
Bowen, G. J. (2018) Gridded maps of the isotopic composition of meteoric waters. http://www.waterisotopes.org.
Bowen, G.J., Wassenaar, L.I. and Hobson, K.A., 2005. Global application of stable hydrogen and oxygen isotopes to wildlife forensics. Oecologia, 143(3), pp.337-348.
}
\examples{
data(d2h_world)
plot(d2h_world)
}

