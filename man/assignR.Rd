\name{assignR}
\alias{assignR}
\title{Tools for inferring geographic origin from isotopic data}
\description{
Routines for rescaling isoscapes using known-origin tissue isotope data, assigning origin of unknown samples, and summarizing and assessing assignment results.
}
\examples{
#### Load library ####
library(assignR)

#### Load North America mask ####
data("naMap")
plot(naMap)

#### Load world precipitation hydrogen isoscape ####
data("d2h_world")
plot(d2h_world)

#### Load hydrogen isotope for human hair in North America ####
d = subOrigData(taxon = "Homo sapiens", reference = "Ehleringer et al. 2008", mask = naMap)

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

#### 26 hypothetical unknown-origin samples ####
id=letters
d2H = seq(-160, -80, by=80/25)
un = data.frame(id,d2H)

#### Assignment for unknown-origin examples ####
asn = pdRaster(r, unknown = un, mask = naMap)

#### SpatialPolygons for two regions of interest ####
data("plover_range_BreedingSeason")
data("plover_range_NonBreedingSeason")
plot(naMap)
lines(plover_range_BreedingSeason, col = c("red"))
lines(plover_range_NonBreedingSeason, col = c("blue"))

#### Get odds ratio for two regions using SpatialPolygon method
p12 = rbind(plover_range_BreedingSeason, plover_range_NonBreedingSeason)
oddsRatio(asn, p12)                      

#### Create SpatialPoints for two points of interest ####
pp1 <- c(-108,42)
pp2 <- c(-103,25)
pp12 <- SpatialPoints(coords = rbind(pp1,pp2), 
  proj4string=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(pp12, pch = 21, bg = "light blue")

#### Get odds ratio for two points using SpatialPoints method ####
oddsRatio(asn, pp12)

#### Binary reclassification ####
# Top 10% of probability surface (defined by % area)
qtlRaster(asn, threshold = 0.1)

# Top 10% of probability surface (defined by % probability)
qtlRaster(asn, threshold = 0.1, thresholdType = "prob")

#### Joint probability for individuals of common origin ####
jointP(asn)

#### Probability that at least one individual came from the location (union of probabilities) ####
unionP(asn)

#### Quality analysis of geographic assignment ####
# get hydrogen isotope values for some known-origin samples
bird_d2h = subOrigData(taxon = "Setophaga ruticilla", mask = naMap)

# run quality assessment using precipitation hydrogen isoscape
d2h_QA <- QA(isoscape = d2h_world, known = bird_d2h, valiStation = 2,
             valiTime = 5, mask = naMap)

# plot the QA result
plot(d2h_QA)
}
