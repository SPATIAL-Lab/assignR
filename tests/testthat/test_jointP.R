context("pdRaster")
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

test_that("pdRaster can correctly calculate posterior probabilities of origin 
          for a sample based on its isotope ratio",{
            expect_equal(cellStats(j, sum), 1)
            expect_is(j, "RasterLayer")
            expect_error(jointP(d))
})