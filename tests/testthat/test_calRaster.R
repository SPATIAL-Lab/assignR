context("calRaster")
library(assignR)
library(raster)
data("naMap")
data("d2h_lrNA")
d = subOrigData(group = "Modern human")
d_hasNA = d
d_hasNA$data$d2H[1] = NA
d_diffProj = d
d_diffProj$data = spTransform(d$data, "+init=epsg:28992")
d_usr_bad = d$data
d_usr_good = d_usr_bad
d_usr_good@data = data.frame(d$data$d2H, d$data$d2H.sd)
d_noCRS = d
crs(d_noCRS$data) = NA

d2h_lrNA_noCRS = d2h_lrNA
crs(d2h_lrNA_noCRS) = NA

mask_diffProj = spTransform(naMap, "+init=epsg:28992")

mask_noCRS = naMap
crs(mask_noCRS) = NA

tempVals = getValues(d2h_lrNA)
tempVals[is.nan(tempVals)] = 9999
d2h_lrNA_with9999 = setValues(d2h_lrNA, tempVals)

d2h_lrNA_na = crop(d2h_lrNA, naMap)

r1 = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap, genplot = FALSE)
r2 = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap, interpMethod = 1, genplot = FALSE)
r3 = calRaster(known = d, isoscape = d2h_lrNA_with9999, NA.value = 9999, genplot = FALSE)

test_that("calRaster can correctly uses known-origin tissue data to rescale a map of 
          environmental isotope values to a map of tissue value (and associated uncertainty) 
          using a linear regression model.",{
        expect_is(r1, "rescale")
        expect_is(r2, "rescale")
        expect_is(r3, "rescale")
        expect_is(calRaster(known = d_usr_good, isoscape = d2h_lrNA), "rescale")
        expect_output(calRaster(known = d, isoscape = d2h_lrNA, outDir = "temp"))
        expect_equal(nlayers(r1$isoscape.rescale), 2)
        expect_error(calRaster(known = d$data$d2H, isoscape = d2h_lrNA))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA, outDir = 2))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA, interpMethod = 3))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA, genplot = 2))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA_noCRS))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA$mean))
        expect_error(calRaster(known = d_usr_bad, isoscape = d2h_lrNA))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA, mask = mask_noCRS))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA, mask = d))
        expect_error(calRaster(known = d_noCRS, isoscape = d2h_lrNA))
        expect_error(calRaster(known = d_hasNA, isoscape = d2h_lrNA, ignore.NA = F))
        expect_error(calRaster(known = d, isoscape = d2h_lrNA_na, ignore.NA = FALSE))
        expect_message(calRaster(known = d_diffProj, isoscape = d2h_lrNA))
        expect_warning(calRaster(known = d, isoscape = d2h_lrNA, mask = mask_diffProj))
        expect_warning(calRaster(known = d, isoscape = d2h_lrNA_na))
})
