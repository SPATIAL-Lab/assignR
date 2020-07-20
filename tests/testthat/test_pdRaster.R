context("pdRaster")
library(assignR)
library(raster)
data("naMap")
data("d2h_world")

d = subOrigData(group = "Modern human", mask = naMap)
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)
id = "smile"
d2H = -80
un = data.frame(id, d2H)
asn = pdRaster(r, unknown = un, mask = naMap)

mask_noCRS = naMap
crs(mask_noCRS) = NA

mask_diffProj = spTransform(naMap, "+init=epsg:28992")

un_hasNA = un
un_hasNA[1,2] = NA

#prior = r$isoscape.rescale$mean
#prior = projectRaster(prior, crs = "+init=epsg:28992")

test_that("pdRaster can correctly calculate posterior probabilities of origin 
          for a sample based on its isotope ratio",{
            expect_is(asn, "RasterLayer")
            expect_error(pdRaster(r$lm.model, unknown = un))
            expect_error(pdRaster(r$isoscape.rescale$mean, unknown = un))
            expect_error(pdRaster(stack(r$isoscape.rescale,r$isoscape.rescale), un))
            expect_error(pdRaster(r, unknown = as.matrix(un)))
            expect_error(pdRaster(r, unknown = un_hasNA))
            expect_error(pdRaster(r, unknown = data.frame(un$id,un$id)))
            expect_error(pdRaster(r, unknown = data.frame(un,un)))
            expect_error(pdRaster(r, unknown = un, genplot = 2))
            expect_error(pdRaster(r, unknown = un, outDir = 2))
            expect_error(pdRaster(r, unknown = un, mask = mask_noCRS))
            expect_error(pdRaster(r, unknown = un, mask = 2))
            expect_warning(pdRaster(r, unknown = un, mask = mask_diffProj))
            
})