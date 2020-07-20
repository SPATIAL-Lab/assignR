context("subOrigData")
library(assignR)
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
d2 = subOrigData(group = c("Passerine"), std_scale = "US_H_5")

test_that("suOrigData can correctly subset known-origin dataset in the assignR package",{
  expect_equal(class(d1), "subOrigData")
  expect_equal(length(d1$data), 188)
  expect_is(d2$data, "SpatialPointsDataFrame")
  expect_error(subOrigData(taxon = "Turdus philomelos", mask = naMap))
  expect_error(subOrigData(taxon = "Turdus philomelos", marker = "d14C"))
  expect_warning(subOrigData(taxon = "Charadrius montanus", age_code = c("chick", "newborn")))
  expect_warning(subOrigData(taxon = c("Danaus plexippus", "Vanellus malabaricus")))
  expect_warning(subOrigData(group = c("Insect","Badgers")))
  expect_warning(subOrigData(reference = c(3, "Ma 2020")))
  expect_warning(subOrigData(reference = c(1, 100)))
})
  