# assignR

## Basic

Data and tools supporting geographic assignment of materials based on their isotopic chemistry. Isoscapes (environmental isotope maps) can be generated externally or defaults (*d2h_world.rda*, *d18o_world.rda*) are provided in the package. Data from samples of known origin are used to calibrate the relationship between isoscape and sample values, and can be provided by the user or extracted from the package database (*knownOrig.rda*). Functions (*calRaster*, *pdRaster*) support calibrating the isoscape and inverting Bayes theorem to estimate the probability of origin for unknown samples across a geographic study domain. Functions (*QA*, *plot.QA*) allow quality assessment of assignment results and comparison of methods using split-sample tests and known origin data. Functions (*oddsRatio*, *qtlRaster*, *jointP*, *unionP*) support post-hoc classification of results, summarization of results from multiple samples, and comparison of support for different locations.

For step-by-step examples demonstrating all functions see [https://spatial-lab.github.io/assignR/](https://spatial-lab.github.io/assignR/).

## Install and load
install.packages("assignR")     
library(assignR)

## Package contents

**Datasets**

*d2h_world.rda* - Global growing season precipitation H isoscape from waterisotopes.org, including predicted mean and 95% confidence interval width

*d2h_lrNA.rda* - Low-resolution, North American crop of d2h_world.rda, used in examples.

*d18o_world.rda* - Global growing season precipitation O isoscape from waterisotopes.org, including predicted mean and 95% confidence interval width

*knownOrig.rda*	- Hydrogen and oxygen isotope values of known-origin samples including human hair, insect chitin and bird feathers, with location information (currently 4218 samples); information on calibration scales used to report data from different labs, useful for converting between scales

*naMap.rda* - North America outline

*states.rda* - 48 contiguous United States

**Functions**

*subOrigData* - Subset the known-origin stable isotope dataset included in this package

*refTrans* - Transform data among calibration scales

*calRaster* - Transform environmental isoscape to tissue isoscape

*pdRaster* - Assign sample to calibrated tissue isoscape based on tissue isotopic composition

*qtlRaster* - Select most likely region of origin from posterior probability surface (by cumulative percent area probability)

*jointP* - Calculate joint probability for individuals of common origin (product of probabilities)

*unionP* - Calculate probability that at least one individual came from each map location (union of probabilities)

*oddsRatio* - Calculate ratio of odds for two locations or areas (points or polygons)

*QA* - Quality analysis of geographic assignment

*plot.QA* - Plot results of one or more quality analyses from QA function

<!-- badges: start -->
  [![Travis build status](https://travis-ci.org/SPATIAL-Lab/assignR.svg?branch=master)](https://travis-ci.org/SPATIAL-Lab/assignR)
  [![codecov](https://codecov.io/gh/Demerara/assignR/branch/master/graph/badge.svg)](https://codecov.io/gh/Demerara/assignR) 
  <!-- badges: end -->

