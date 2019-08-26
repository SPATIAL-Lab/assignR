\name{oddsRatio}
\alias{oddsRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Odds ratio of points or regions
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

#### Load hydrogen isotope for human hair in North America ####
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap)

# Create SpatialPolygons for two regions of interest
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

# Create SpatialPoints with two points of interest
pp1 <- c(-100,45)
pp2 <- c(-100,60)
pp12 <- SpatialPoints(coords = rbind(pp1,pp2))

# Caculate odds ratio for the two points created above
oddsRatio(asn, pp12)
}

