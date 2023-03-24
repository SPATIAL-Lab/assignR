library(testthat)
library(assignR)
library(terra)
library(sp)

d2h_lrNA = rast(system.file("extdata/d2h_lrNA.tif", package = "assignR"))
sr_MI = rast(system.file("extdata/sr_MI.tif", package = "assignR"))

test_check("assignR")
