\name{oddsRatio}

\alias{oddsRatio}

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
SpatRaster of probability density maps, e.g., as produced by \code{pdRaster}.
}
  \item{inputP}{
SpatialPoints* object of length 1 or 2 or SpatialPolygons* object of length 2
}
}

\examples{
library(sp)
library(terra)

# load North America boundary and global isoscape
data("naMap")
d2h_lrNA = rast(system.file("extdata/d2h_lrNA.tif", package = "assignR"))

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, genplot = FALSE)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(d, d2h_lrNA, naMap, genplot = FALSE)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id, d2H)

# posterior probabilities
pp = pdRaster(r, un, mask = naMap, genplot = FALSE)

# SpatialPolygons for two regions of interest
data("states")
s1 = states[states$STATE_ABBR == "UT",]
s2 = states[states$STATE_ABBR == "NM",]
plot(naMap)
plot(s1, border = "red", add = TRUE)
plot(s2, border = "blue", add = TRUE)

# Get odds ratio for two regions using SpatialPolygon method
s12 = rbind(s1, s2)
oddsRatio(pp, s12)                     

# Create SpatialPoints for two points of interest
p1 = c(-112, 40)
p2 = c(-105, 33)
p12 = SpatialPoints(coords = rbind(p1, p2), 
  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
points(p12, pch = 21, bg = "light blue")

# Get odds ratio for two points using SpatialPoints method
oddsRatio(pp, p12)
}

