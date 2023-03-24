iso = getIsoscapes("CaribSr")
iso2 = getIsoscapes("CaribSr")

test_that("getIsoscapes works",{
<<<<<<< HEAD
  if(!is.null(iso)){
    expect_is(iso, "SpatRaster")
  }
=======
  expect_is(iso, "SpatRaster")
>>>>>>> c37cfc69cc00426bb1adaa411491964227e16adb
  expect_equal(iso, iso2)
  expect_error(getIsoscapes("USPrecipALLs"))
})
  