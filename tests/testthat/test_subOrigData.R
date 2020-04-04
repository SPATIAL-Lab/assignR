context("subOrigData")
library(assignR)
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
d2 = subOrigData(group = c("insect","bird"))

test_that("suOrigData can correctly subset known-origin dataset in the assignR package",{
  expect_equal(length(d1), 133)
  expect_is(d2, "SpatialPointsDataFrame")
  expect_error(subOrigData(taxon = "Turdus philomelos", mask = naMap))
  expect_error(subOrigData(taxon = "Turdus philomelos", marker = "d14C"))
  expect_warning(subOrigData(taxon = "Charadrius montanus", age_code = c("chick", "newborn")))
  expect_warning(subOrigData(taxon = c("Danaus plexippus", "Vanellus malabaricus")))
  expect_warning(subOrigData(group = c("insect","Badgers")))
  expect_warning(subOrigData(reference = c("Wunder 2007", "Ma 2020")))
})
  