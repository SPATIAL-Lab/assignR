context("isoStack")
library(assignR)
library(raster)
data("d2h_lrNA")

d2h_crop = crop(d2h_lrNA, extent(matrix(c(-100, -60, 25, 45), nrow = 2)))
d2h_aea = projectRaster(d2h_lrNA, crs = crs("+proj=aea +lat_1=20 
                                            +lat_2=60 +lat_0=40 
                                            +lon_0=-96 +x_0=0 +y_0=0 
                                            +ellps=GRS80 +datum=NAD83"))

test_that("isoStack can stack isoscapes",{
  expect_error(isoStack(d2h_lrNA[[1]], d2h_lrNA))
  expect_error(isoStack(d2h_lrNA))
  expect_error(isoStack(d2h_lrNA, d18o_world, clean = FALSE))
  expect_s3_class(isoStack(d2h_lrNA, d2h_crop), "isoStack")
  expect_s3_class(isoStack(d2h_lrNA, d2h_aea), "isoStack")
})