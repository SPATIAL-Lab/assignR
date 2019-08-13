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
qtlRaster(pdR, threshold, thresholdType)
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
numeric 1 or 2. 1 is to choose top percentage by cumulative probability; 2 is to choose top percentage by area
}
  \item{genplot}{
generate summary plot? (T or F)
}
  \item{pdf}{
save summary plot as PDF? (T or F)
}
}

\examples{
# load North America SpatialPolygonsDataFrame data in assignR package
data("naMap")
data(d2h_world)
d = subOrigData(taxon = c("Homo sapiens"), mask = naMap)
# Exclude some outliers. This step is optional, which depends on your data quality
d <-as.data.frame(d)
dd = d[d$coords.x1<(-80),]
dd <- SpatialPointsDataFrame(dd[,2:3], as.data.frame(dd[,1]))
crs(dd) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# Rescale from environmental isoscape to tissue isoscape (known must be SpatialPointsDataFrame with coordinate reference system )
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)
# four unknown-origin example
id = c("A", "B", "C", "D")
d2H = c(-110, -180, -130, -150)
un = data.frame(id,d2H)
# assignment
asn = pdRaster(r,unknown=un,mask=naMap)
# top 10% of probability surface (defined by % area)
qtlRaster(asn, threshold = 0.1, thresholdType = 2, pdf = T)
# top 10% of probability surface (defined by % cumulative probability)
qtlRaster(asn, threshold = 0.1, thresholdType = 1, pdf = F)
}
