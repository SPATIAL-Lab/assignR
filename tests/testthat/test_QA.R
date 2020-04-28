context("QA")
library(assignR)
library(raster)
data("naMap")
data("d2h_world")
data("knownOrig")
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")
qa1 = QA(isoscape = raster::aggregate(d2h_world, 6), known = d1, 
         valiStation = 1, valiTime = 2, mask = naMap, name = "Charadrius")
qa2 = QA(isoscape = raster::aggregate(d2h_world, 6), known = d2, 
         valiStation = 1, valiTime = 2, mask = naMap, name = "Buteo", setSeed = F)

# d2h_world_noCRS = d2h_world
# crs(d2h_world_noCRS) = NA
# 
# d2_hasNA = d2
# d2_hasNA[1,1] = NA
# 
# d2_noCRS = d2
# crs(d2_noCRS) = NA
# 
# d2_diffProj = spTransform(d2, "+init=epsg:28992")
# 
# d2_2columData = d2
# d2_2columData@data$col2 = d2_2columData@data$d2H
# 
# mask_noCRS = naMap
# crs(mask_noCRS) = NA
# 
# mask_diffProj = spTransform(naMap, "+init=epsg:28992")

test_that("QA can correctly evaluate how well does a given isoscape and known origin data set 
          constrain the geographic origin of samples",{
            expect_is(qa1, "QA")
            expect_is(qa2, "QA")
            
            #expect_error(QA(d2h_world_noCRS, known = d2, valiStation = 1))
            #expect_error(QA(d2h_world$mean, known = d2, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2@data, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2, name = 1, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2_hasNA, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2_noCRS, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2_2columData, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2, valiStation = 100))
            #expect_error(QA(d2h_world, known = d2, valiTime = 1, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2, mask = 2, valiStation = 1))
            #expect_error(QA(d2h_world, known = d2, mask = mask_noCRS, valiStation = 1))
            
            #expect_warning(QA(d2h_world, known = d2, valiStation = 1, mask = mask_diffProj))
            #expect_warning(QA(d2h_world, known = d2_diffProj, valiStation = 1, mask = naMap))
   
})