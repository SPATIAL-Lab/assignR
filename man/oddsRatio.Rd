\name{oddsRatio}
\alias{oddsRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
oddsRatio
}
\description{
Calculate ratio of odds for two locations (points or polygons)
}
\usage{
oddsRatio(pdR, inputP)
}
\arguments{
  \item{pdR}{
probability density map which must be raster or RasterLayer
}
  \item{inputP}{
points or polygons which must be data.frame or SpatialPolygons
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

# Create SpatialPolygons with two polygons
p1 <- c(-100,60,-100,65,-110,65,-110,60,-100,60)
p1 <- matrix(p1, 5,2, byrow = T)
p1 <- Polygon(p1)
p1 <- Polygons(list(p1), "p1")
p2 <- c(-100,40,-100,45,-110,45,-110,40,-100,40)
p2 <- matrix(p2, 5,2, byrow = T)
p2 <- Polygon(p2)
p2 <- Polygons(list(p2), "p2")
p12 <- SpatialPolygons(list(p1,p2),1:2)
plot(p12)
# Calculate odds ratio for the two polygons created above
oddsRatio(asn, p12)

# points
pp1 <- c(-100,45)
pp2 <- c(-100,60)
pp12 <- as.data.frame(rbind(pp1,pp2))
# Caculate odds ratio for the two points created above
oddsRatio(asn, pp12)
}

