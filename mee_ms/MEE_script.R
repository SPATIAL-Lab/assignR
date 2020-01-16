#Load packages
library(assignR)

#Set up data and projection
d = subOrigData(taxon = "Lanius ludovicianus", reference = "Hobson et al. 2012", mask = naMap)
s = readOGR("states_shapefile/states.shp")
p = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83"

#Project our spatial objects
isona = crop(d2h_world, naMap)
isona = projectRaster(isona, crs = p)
s = spTransform(s, p)
c = spTransform(naMap, p)

#Calibrate tissue isoscape
r = calRaster(known = d, isoscape = isona, mask = naMap)

#Hypothetical unknown-origin sample
un = data.frame("id" = "A", "d2H" = -110)

#Posterior probability raster
pdR = pdRaster(r, unknown = un)

#Two states to use in comparing odds
s1 = s[s$STATE_ABBR == "UT",]
s2 = s[s$STATE_ABBR == "NM",]
s12 = rbind(s1, s2)

#Calculate odds ratio
oddsRatio(pdR, s12)

#Upper 90% probability region
q = qtlRaster(pdR, 0.9, "prob")

#Map figure
png("Figure2.png", width = 10, height = 7.5, units = "in", res = 600)
layout(matrix(c(1,2), nrow = 1))
plot(pdR*1000, xlim = c(-4.5e6, 3.5e6), ylim = c(-2e6, 5e6), axes = FALSE, box = FALSE)
lines(c)
lines(s1, col = "red")
lines(s2, col = "blue")
text(-4.2e6, 5e6, "(a)")
text(3.3e6, 1.2e6, "Probability (x 1000)", srt = 90)

plot(q, xlim = c(-4.5e6, 3.5e6), ylim = c(-2e6, 5e6), axes = FALSE, box = FALSE, legend = FALSE)
lines(c)
text(-4.2e6, 5e6, "(b)")

dev.off()

#Generate random points 
locs = spsample(c, 300, "random")

#Get isotope values
locs$vals = extract(isona[[1]], locs)
locs = locs[!is.na(locs$vals),]
locs = locs[1:200,]

#Add some reasonable random noise to the known origin data
locs$vals = locs$vals + rnorm(nrow(locs), 0, 10)
k = SpatialPointsDataFrame(locs, data.frame(d2H = locs$vals), proj4string = p)

#Create a version of the isoscape that has no uncertainty
isoperf = isona
isoperf[[2]] = setValues(isoperf[[2]], 0)

#Now assign to the perfect raster; this is 'correct' since the isoscape values
#used to generate the synthetic data were known perfectly
qa_perf = QA(isoperf, k, name = "Not inflated")

#Now assign using the original isoscape, with error
qa_err = QA(isona, k, name = "Inflated")

#Plot it
plot(qa_perf, qa_err, outDir = "C:/Users/u0133977/Dropbox/Hypomirror/assignR_MEE")

