# assignR news

## assignR 2.2.3.9000
* 

## assignR 2.2.3
* Bug fixes

## assignR 2.2.2
* Updated error handling in getIsoscapes

## assignR 2.2.1
* Remove data from projects 14 and 15 from knownOrig dataset

## assignR 2.2.0
* Add wDist function and c and plot methods for summarizing weighted distance and bearing distributions using sample collection locations and posterior probability maps
* QA option to run without iterative recalibration of isoscape
* Bug fixes
* Documentation edits

## assignR 2.1.1
* Bug fixes

## assignR 2.1.0
* Add isoStack function to stack multiple isoscapes in a single data object; added plot s3 method for this class
* Add getIsoscape function supporting download of gridded isotope maps; removed global precipitation maps previously distributed with package
* Update pdRaster and QA to support multivariate analysis
* Add option to include spatial prior in QA
* Bug fixes

## assignR 2.0.0

* knownOrig database has been expanded and reformatted
* New data objects document different calibration standards used to generate known-origin H and O isotope data
* subOrigData supports transformation of data among different calibration standard scales using the 'standard-chain' method of Magozzi et al. (in press); format of return object from this function has changed
* calRaster changes including new format for input object "known" and use of weighted least squares regression; tissue isoscape variance calculation updated to a + b - c, where a is isoscape grid cell variance, b is the residual variance of tissue predictions made from isoscape-tissue rescaling functions fit using values sampled from the isoscape with noise, and c is the variance of the sampled isoscape values
* QA changes including new format for input object "known" and option to resample known data by site rather than by sample; for bySite option returned results are the average of site-level average statistics; argument order changed for consistency with calRaster
* sp objects updated to support WKT2_2019 strings

## assignR 1.2.1

* Updates to ensure compatibility with new CRS specifications in rgdal
* QA now accepts optional argument "by", allowing reduced run-time

## assignR 1.2.0

* NAMESPACE imports required functions only
* Implemented testthat testing
* Minor corrections to knownOrig data
* Minor enhancements and bug fixes

## assignR 1.1.3.1

* Data update: remove plover maps, add US states

## assignR 1.1.3

* Shorter run-time for vignette examples

## assignR 1.1.2

* Functional examples for QA
* Improved handling of graphical parameters in plotting functions
* Bug fixes

## assignR 1.1.1

* Consistent syntax for saving output to disk, which requires explicit specification of directory
* Improved messaging from functions
* Bug fixes

## assignR 1.1

* Propagation of error covariance in isoscape and rescaling models removes source of bias in posterior probabilities
* plot.QA converted to S3 method
* Bug and code fixes to pass CRAN check

## assignR 1.0

Initial GitHUB release
