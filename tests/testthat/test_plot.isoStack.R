context("plot.isoStack")
library(assignR)
library(raster)
data("d2h_lrNA")

d2h_crop = crop(d2h_lrNA, extent(matrix(c(-100, -60, 25, 45), nrow = 2)))
s = isoStack(d2h_lrNA, d2h_crop)

test_that("plot can plot isoStack",{
  expect_error(plot.isoStack(d2h_lrNA))
  expect_silent(plot(s))
})