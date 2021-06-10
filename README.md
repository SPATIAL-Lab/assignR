# assignR

## Basic

Data and tools supporting geographic assignment of materials based on their isotopic chemistry. Isoscapes (environmental isotope maps) can be generated externally or downloaded using the *getIsoscapes* function. Data from samples of known origin are used to calibrate the relationship between isoscape and sample values, and can be provided by the user or extracted from the package database (*knownOrig.rda*). Database data or user-provided known-origin or unknown origin sample data can be transformed among different H and O isotope reference scales to improve comparability (*refTrans*). Functions (*calRaster*, *pdRaster*) support calibrating one or more isoscapes (multiple layers combined using *isoStack*) and inverting the assignment model to estimate the probability of origin for unknown samples across a geographic study domain. Functions (*QA*, *plot.QA*) allow quality assessment of assignment results and comparison of methods using split-sample tests and known origin data. Functions (*oddsRatio*, *qtlRaster*, *jointP*, *unionP*) support post-hoc classification of results, summarization of results from multiple samples, and comparison of support for different locations.

For current production release, see the vignette [here](https://CRAN.R-project.org/package=assignR) and install from CRAN.

For examples demonstrating functions in the latest development release, see [https://spatial-lab.github.io/assignR/](https://spatial-lab.github.io/assignR/).

## Install and load latest CRAN release
install.packages("assignR")     
library(assignR)

## Package contents

**Datasets**

*d2h_lrNA.rda* - Low-resolution, North American crop of growing season d2H isoscape, used in examples.

*sr_MI.rda* - Low-resolution crop of locally-weathered Sr isoscape, used in examples.

*knownOrig.rda*	- Hydrogen and oxygen isotope values of known-origin samples including human hair, insect chitin and bird feathers, with location information (currently 4218 samples)

*stds.rda* - Information on reference scales used to report data from different labs, useful for converting between scales

*naMap.rda* - North America outline

*states.rda* - 48 contiguous United States

**Functions**

*subOrigData* - Subset the known-origin stable isotope dataset included in this package

*refTrans* - Transform data among reference scales

*getIsoscapes* - Download and unpack isoscapes from waterisotopes.org

*isoStack* - Combine multiple isoscapes

*plot.isoStack* - Plot isoStack object

*calRaster* - Transform one or more isoscapes to reflect target sample type

*pdRaster* - Assign sample to calibrated isoscape(S) based on isotopic composition(s)

*qtlRaster* - Select most likely region of origin from posterior probability surface (by cumulative percent area probability)

*jointP* - Calculate joint probability for individuals of common origin (product of probabilities)

*unionP* - Calculate probability that at least one individual came from each map location (union of probabilities)

*oddsRatio* - Calculate ratio of odds for two locations or areas (points or polygons)

*QA* - Quality analysis of geographic assignment

*plot.QA* - Plot results of one or more quality analyses from QA function

<!-- badges: start -->
  [![Build status](https://github.com/SPATIAL-Lab/assignR/actions/workflows/r.yml/badge.svg)](https://github.com/SPATIAL-Lab/assignR/actions)
  [![codecov](https://codecov.io/gh/SPATIAL-Lab/assignR/branch/master/graph/badge.svg)](https://codecov.io/gh/SPATIAL-Lab/assignR) 
  <!-- badges: end -->

