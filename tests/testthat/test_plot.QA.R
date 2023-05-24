library(terra)
dmp = capture_output({
  r1 = aggregate(sr_MI, 5)
  r2 = project(d2h_lrNA, r1)
  r3 = isoStack(r1, r2)
  dpts = spatSample(c(r1, r2), 20, na.rm = TRUE, as.points = TRUE)
  dpts$weathered.mean = dpts$weathered.mean + rnorm(20, 0, 0.0005)
  dpts$d2h = dpts$d2h + rnorm(20, 0, 8)
  dpts$Site_ID = letters[1:20]
  d1 = dpts[,1:2]
  d2 = dpts
  d3 = suppressWarnings(subOrigData(group = c("Raptor", "Passerine", "Water bird"), 
                                    mask = states[states$STATE_ABBR == "MI",],
                                    ref_scale = NULL, genplot = FALSE))
  
  qa1 = suppressWarnings(QA(d1, r1, valiStation = 1, valiTime = 2, by = 25, 
                            name = "Sr", bySite = FALSE))
  qa2 = suppressWarnings(QA(d2, r3, valiStation = 1, 
                            valiTime = 2, by = 25, 
                            name = "Multi", setSeed = F))
  qa3 = suppressWarnings(QA(d3, r2, valiStation = 1, 
                            valiTime = 2, by = 25, 
                            name = "SOD", setSeed = F))
})

test_that("QA and plot.QA work",{
  expect_is(qa1, "QA")
  expect_is(qa2, "QA")
  expect_silent(plot(qa1, qa2))
  expect_silent(plot(qa3))
  expect_silent(plot(qa1, qa2))
  expect_silent(plot(qa1, outDir = tempdir()))
})
