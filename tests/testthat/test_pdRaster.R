context("pdRaster")
library(assignR)
library(raster)
data("naMap")
data("d2h_lrNA")
data("sr_MI")

d = subOrigData(group = "Modern human", mask = naMap)
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap, genplot = FALSE)
id = "smile"
d2H = -80
d2H.sd = 2
un = data.frame(id, d2H, d2H.sd, d2H_cal = "UT_H_1")

r2 = r
r2[[1]][[1]] = r2[[1]][[1]] + rnorm(ncell(r2[[1]][[1]]), 0, 10)
rmulti = isoStack(r, r2)
un2 = refTrans(un)

r3 = isoStack(r, sr_MI)
un3 = data.frame(id, d2H, "Sr" = 0.710)

mask_noCRS = naMap
crs(mask_noCRS) = NA

mask_diffProj = spTransform(naMap, "+init=epsg:28992")

un_hasNA = un
un_hasNA[1,2] = NA

test_that("pdRaster can correctly calculate posterior probabilities of origin 
          for a sample based on its isotope ratio",{
            expect_warning(pdRaster(r, unknown = un, mask = naMap, genplot = FALSE))
            expect_error(pdRaster(rmulti, list(un, un), mask = naMap, 
                                         genplot = FALSE))
            expect_is(pdRaster(rmulti, list(un2, un2), mask = naMap, 
                               genplot = FALSE), "RasterLayer")
            expect_is(pdRaster(r3, un3), "RasterLayer")
            expect_error(pdRaster(r$lm.model, unknown = un))
            expect_error(pdRaster(r$isoscape.rescale$mean, unknown = un))
            expect_error(pdRaster(stack(r$isoscape.rescale,r$isoscape.rescale), un))
            expect_error(pdRaster(r, unknown = as.matrix(un)))
            expect_error(pdRaster(r, unknown = un_hasNA))
            expect_error(pdRaster(r, unknown = data.frame(un$id,un$id)))
            expect_error(pdRaster(r, unknown = data.frame(un[,1])))
            expect_error(pdRaster(r, unknown = un, genplot = 2))
            expect_error(pdRaster(r, unknown = un, outDir = 2))
            expect_error(pdRaster(r, unknown = un, mask = mask_noCRS))
            expect_error(pdRaster(r, unknown = un, mask = 2))
            expect_message(pdRaster(r, unknown = un, mask = mask_diffProj, genplot = FALSE))
})