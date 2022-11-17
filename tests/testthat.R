library(testthat)
library(assignR)
library(raster)
library(terra)

d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))
sr_MI = rast(system.file("data/sr_MI.tif", package = "assignR"))

test_check("assignR")