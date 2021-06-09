library(assignR)
data("naMap")
data("d2h_lrNA")
d = subOrigData(taxon = "Homo sapiens", 
                dataset = 10, mask = naMap)
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)
asn = pdRaster(r, unknown = un, mask = naMap)
j = jointP(asn)
u = unionP(asn)

data("states")
s1 = states[states$STATE_ABBR == "UT",]
s2 = states[states$STATE_ABBR == "NM",]
s12 = rbind(s1, s2)
o1 = oddsRatio(asn, s12)                     
pp1 = c(-112,40)
pp2 = c(-105,33)
pp12 = SpatialPoints(coords = rbind(pp1,pp2), 
                     proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
o2 = oddsRatio(asn, pp12)
o3 = oddsRatio(asn, pp12[1])
o4 = oddsRatio(asn$A, pp12)
o5 = oddsRatio(asn$A, s12)

s12_diffProj = spTransform(s12, CRS("+init=epsg:28992"))
pp12_diffProj = spTransform(pp12, CRS("+init=epsg:28992"))

pp12_noCRS = pp12
crs(pp12_noCRS) = NA

s12_noCRS = s12
crs(s12_noCRS) = NA

pp3 = pp1
pp121 = SpatialPoints(coords = rbind(pp1,pp2,pp3), 
                      proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

test_that("post-processing tools work",{
            expect_equal(cellStats(j, sum), 1)
            expect_is(j, "RasterLayer")
            expect_error(jointP(d))
            
            expect_is(u, "RasterLayer")
            expect_error(unionP(d2H))
            
            expect_is(o1, "list")
            expect_is(o2, "list")
            expect_is(o3, "data.frame")
            expect_is(o4, "list")
            expect_is(o5, "list")
            
            expect_error(oddsRatio(naMap,s12))
            expect_error(oddsRatio(asn, data.frame(30.6, 50.5)))
            expect_error(oddsRatio(asn, s12_noCRS))
            expect_error(oddsRatio(asn, pp121))
            expect_error(oddsRatio(asn, s1))
            expect_error(oddsRatio(asn, pp12_noCRS))
            
            expect_message(oddsRatio(asn, s12_diffProj))
            expect_message(oddsRatio(asn, pp12_diffProj))
})