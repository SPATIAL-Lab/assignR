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

#### Exclude some outliers. This step is optional, which depends on your data quality ####
d <-as.data.frame(d)
dd = d[d$coords.x1<(-80),]
dd <- SpatialPointsDataFrame(dd[,2:3], as.data.frame(dd[,1]))
crs(dd) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#### Rescale from environmental isoscape to tissue isoscape (known must be SpatialPointsDataFrame with coordinate reference system ) ####
r = calRaster(known = dd, isoscape = d2h_world, mask = naMap)

#### 26 unknown-origin examples ####
id=letters
d2H = seq(-160, -80, by=80/25)
un = data.frame(id,d2H)

#### Assignment for unknown-origin examples ####
asn = pdRaster(r,unknown=un,mask=naMap)

#### Create SpatialPolygons with two polygons ####
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

#### Create data.frame with two points ####
pp1 <- c(-100,45)
pp2 <- c(-100,60)
pp12 <- as.data.frame(rbind(pp1,pp2))
#### Caculate odds ratio for the two polygons created above ####
oddsRatio(asn, p12)

#### Caculate odds ratio for the two points created above ####
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
# oxygen and hydrogen isotopes of known-origin bird
data(bird_isotope) 

# crop the world hydrogen data to North America
r <- crop(d2h_world, naMap)
plot(r)

# convert 2 standard deviation from d2h_world to 1 standard deviation
r[[2]] <- r[[2]]/2

# seperate the hydrogen isotope for the known-origin bird
bird_d2h <- bird_isotope[1:20,c("Longitude", "Latitude", "d2H")]
coordinates(bird_d2h) <- c(1,2)
proj4string(bird_d2h) <- proj4string(d2h_world)

# run quality assessment based hydrogen isotope from precipitation and known-origin bird
d2h_QA <- QA(isoscape = r, known = bird_d2h, valiStation = 2,
             valiTime = 5, setSeed = T)

# plot the QA result
plot(d2h_QA)
}
