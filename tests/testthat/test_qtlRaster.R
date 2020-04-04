context("qtlRaster")
library(assignR)
data("naMap")
data("d2h_world")
d = subOrigData(taxon = "Homo sapiens", 
                reference = "Ehleringer et al. 2008", mask = naMap)
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id, d2H)
asn = pdRaster(r,unknown=un,mask=naMap)
q1 = qtlRaster(asn, threshold = 0.1, thresholdType = "area")
q2 = qtlRaster(asn, threshold = 0.1, thresholdType = "prob")
q3 = qtlRaster(asn, threshold = 0)

test_that("qtlRaster can correctly: Selects the grid cells of 
          probability density rasters with the highest probability 
          and returns rasters with these cell values set to 1;
          Cells are selected based on the user-specified quantile threshold 
          so that the most-probable cells representing a given fraction 
          of the assignment area or posterior probability are returned.",{
            expect_is(q1, "RasterStack")
            expect_equal(nlayers(q1), 4)
            expect_equal(nlayers(q2), 4)
            expect_equal(nlayers(q3), 4)
            expect_error(qtlRaster(asn, threshold = "a"))
            expect_error(qtlRaster(asn, threshold = 10))
            expect_error(qtlRaster(asn, threshold = "a"), thresholdType = "probability")
            expect_error(qtlRaster(asn, threshold = 0.1, genplot = "A"))
            expect_error(qtlRaster(asn, threshold = 0.1, outDir = 1))
            
})
