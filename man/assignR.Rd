\name{assignR}
\alias{assignR}

\title{assignR}

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
d = subOrigData(taxon = c("Homo sapiens"), mask = naMap)

#### Exclude some outliers. This step is optional, which depends on your data quality.
dd = d[d@coords[,1] < -80,]

#### Rescale from environmental isoscape to tissue isoscape ####
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)

#### 26 hypothetical unknown-origin samples ####
id=letters
d2H = seq(-160, -80, by=80/25)
un = data.frame(id,d2H)

#### Assignment for unknown-origin examples ####
asn = pdRaster(r,unknown=un,mask=naMap, saveFile=FALSE)

#### Create SpatialPolygons for two regions of interest ####
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

#### Create data.frame with two points of interest ####
pp1 <- c(-100,45)
pp2 <- c(-100,60)
pp12 <- as.data.frame(rbind(pp1,pp2))

#### Caculate odds ratio for each sample for the two polygons created above ####
oddsRatio(asn, p12)

#### Caculate odds ratio for each sample for the two points created above ####
oddsRatio(asn, pp12)

#### Binary reclassification ####
### Top 10% of probability surface (defined by % area)
qtlRaster(asn, threshold = 0.1, thresholdType = 2)
### Top 10% of probability surface (defined by % cumulative probability)
qtlRaster(asn, threshold = 0.1, thresholdType = 1)

#### Joint probability for individuals of common origin ####
jointP(asn)

#### Probability that at least one individual came from the location (union of probabilities) ####
unionP(asn)

#### Quality analysis of geographic assignment ####
# get hydrogen isotope values for some known-origin samples
bird_d2h = subOrigData(taxon = "Setophaga ruticilla", mask = naMap)

# run quality assessment using precipitation hydrogen isoscape
d2h_QA <- QA(isoscape = d2h_world, known = bird_d2h, valiStation = 2,
             valiTime = 5, mask = naMap, setSeed = TRUE)

# plot the QA result
plot(d2h_QA)
}
