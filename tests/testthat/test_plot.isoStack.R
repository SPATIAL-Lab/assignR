d2h_crop = crop(d2h_lrNA, ext(-100, -60, 25, 45))
d2h_proj = project(d2h_lrNA, "+proj=longlat +ellps=clrk66 
                   +datum=NAD27 +no_defs")
s = suppressWarnings(isoStack(d2h_lrNA, d2h_crop))

test_that("isoStack can stack isoscapes and plot can plot them",{
  expect_error(isoStack(d2h_lrNA[[1]], d2h_lrNA))
  expect_error(isoStack(d2h_lrNA))
  expect_error(isoStack(d2h_lrNA, d2h_proj, clean = FALSE))
  expect_s3_class(s, "isoStack")
  expect_s3_class(suppressWarnings(isoStack(d2h_lrNA, d2h_proj)), "isoStack")
  expect_error(plot.isoStack(d2h_lrNA))
  expect_silent(plot(s))
})
