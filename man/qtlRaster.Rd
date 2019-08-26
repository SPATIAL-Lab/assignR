\name{qtlRaster}
\alias{qtlRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary reclassification
}
\description{
Select top X percent of probability surface (defined by area percentage or cumulative probability)
}
\usage{
qtlRaster(pdR, threshold, thresholdType = 2, genplot = TRUE, savePDF = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdR}{
probability density map which must be raster or RasterLayer
}
  \item{threshold}{
numeric from 0 to 1, top percentage to choose
}
  \item{thresholdType}{
numeric 1 or 2. 1 selects the top percentage by probability; 2 selects the top percentage by area
}
  \item{genplot}{
generate summary plot? (TRUE or FALSE)
}
  \item{savePDF}{
save summary plot as PDF? (TRUE or FALSE)
}
}

\examples{
# load North America SpatialPolygonsDataFrame data in assignR package
data("naMap")
data(d2h_world)

#### Load hydrogen isotope for human hair in North America ####
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r,unknown=un,mask=naMap)

# top 10% of probability surface (defined by % area)
qtlRaster(asn, threshold = 0.1, savePDF = TRUE)

# top 10% of probability surface (defined by % probability)
qtlRaster(asn, threshold = 0.1, thresholdType = 1)
}
