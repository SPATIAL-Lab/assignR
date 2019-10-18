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
RasterStack or RasterBrick of probability density maps, e.g., as produced by \code{calRaster}. All layers must have common projection.
}
  \item{inputP}{
SpatialPoints or SpatialPolygons (or *DataFrame equivalent) of length 2
}
}

\examples{
# load North America boundary and global isoscape
data("naMap")
data("d2h_world")

# load hydrogen isotope data for human hair in North America
d = subOrigData(taxon = "Homo sapiens", 
  reference = "Ehleringer et al. 2008", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# assignment
asn = pdRaster(r, unknown = un, mask = naMap)

# SpatialPolygons for two regions of interest
data("plover_range_BreedingSeason")
data("plover_range_NonBreedingSeason")
plot(naMap)
plot(plover_range_BreedingSeason, col = c("red"), add = TRUE)
plot(plover_range_NonBreedingSeason, col = c("blue"), add = TRUE)

# Get odds ratio for two regions using SpatialPolygon method
p34 = rbind(plover_range_BreedingSeason, plover_range_NonBreedingSeason)
oddsRatio(asn, p34)                      

# Create SpatialPoints for two points of interest
pp1 = c(-108,42)
pp2 = c(-103,25)
pp12 = SpatialPoints(coords = rbind(pp1,pp2), 
  proj4string=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(pp12, pch = 21, bg = "light blue")

# Get odds ratio for two points using SpatialPoints method
oddsRatio(asn, pp12)
}

