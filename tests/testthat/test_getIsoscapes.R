iso = getIsoscapes("CaribSr")
iso2 = getIsoscapes("CaribSr")

test_that("getIsoscapes works",{
  expect_is(iso, "SpatRaster")
  expect_equal(iso, iso2)
  expect_error(getIsoscapes("USPrecipALLs"))
})
  