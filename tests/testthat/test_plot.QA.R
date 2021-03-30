context("plot.QA")
library(assignR)
data("naMap")
data("d2h_lrNA") 
data("knownOrig") 
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")
qa1 = QA(known = d1, isoscape = d2h_lrNA, valiStation = 1, valiTime = 2, 
         by = 25, mask = naMap, name = "Charadrius")
qa2 = QA(known = d2, isoscape = d2h_lrNA, valiStation = 1, valiTime = 2, 
         by = 25, mask = naMap, name = "Buteo")

test_that("plot.QA can correctly plot the output from QA",{
            expect_error(plot.QA(d1))
            expect_error(plot.QA(qa1, outDir = 2))
            expect_silent(plot.QA(qa1, qa2))
            expect_silent(plot.QA(qa1))
})