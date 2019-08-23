\name{d18o_world}
\alias{d18o_world}
\docType{data}
\title{
data: d18o_world
}
\description{
Global growing season precipitation δ18O isoscape from waterisotopes.org
}
\usage{data("d18o_world")}
\format{
 raster with two layers. The first layer is the mean prediction and the second one is 1 standard deviation
}
\references{
Bowen, G. J. (2018) Gridded maps of the isotopic composition of meteoric waters. http://www.waterisotopes.org.
Bowen, G.J., Wassenaar, L.I. and Hobson, K.A., 2005. Global application of stable hydrogen and oxygen isotopes to wildlife forensics. Oecologia, 143(3), pp.337-348.
}
\examples{
data(d18o_world)
plot(d18o_world)
}
