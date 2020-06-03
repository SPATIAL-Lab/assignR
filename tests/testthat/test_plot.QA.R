context("plot.QA")
library(assignR)
data("naMap")
data("d2h_world") 
data("knownOrig") 
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")
qa1 = QA(isoscape = d2h_lrNA, known = d1, valiStation = 1, valiTime = 2, 
         by = 25, mask = naMap, name = "Charadrius")
plot(qa1)
qa2 = QA(isoscape = d2h_lrNA, known = d2, valiStation = 1, valiTime = 2, 
         by = 25, mask = naMap, name = "Buteo")
plot(qa1, qa2)

test_that("plot.QA can correctly plot the output from QA",{
            expect_error(plot.QA(d1))
            expect_error(plot.QA(qa1, outDir = 2))
            expect_silent(plot.QA(qa1, outDir = "temp"))
            expect_silent(plot.QA(qa1, qa2, outDir = "temp"))
            
})