context("QA")
library(assignR)
library(raster)
data("naMap")
data("d2h_world")
data("knownOrig")
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")
qa1 = QA(isoscape = d2h_lrNA, known = d1, valiStation = 1, valiTime = 2, by = 25, 
         mask = naMap, name = "Charadrius")
qa2 = QA(isoscape = d2h_lrNA, known = d2, valiStation = 1, valiTime = 2, by = 25, 
         mask = naMap, name = "Buteo", setSeed = F)

test_that("QA can correctly evaluate how well a given isoscape and known origin 
          data set constrains the geographic origin of samples",{
            expect_is(qa1, "QA")
            expect_is(qa2, "QA")
})