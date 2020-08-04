context("unionP")
library(assignR)
data("naMap")
data("d2h_world")
d = subOrigData(taxon = "Homo sapiens", dataset = 10, mask = naMap)
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)
asn = pdRaster(r, unknown = un, mask = naMap)
u = unionP(asn)

test_that("unionP can correctly calculate probabilities 
          that at least one individual came from each location 
          in the assignment area (union of probabilities)",{
            expect_is(u, "RasterLayer")
            expect_error(unionP(d2H))
})