# assignR news

## assignR 1.2.1.9000

* knownOrig database has been expanded and reformatted
* New data objects document different calibration standards used to generate known-origin H and O isotope data
* subOrigData supports transformation of data among different calibration standard scales using the 'standard-chain' method of Magozzi et al. (in prep); format of return object from this function has changed
* calRaster changes including new format for input object "known" and use of weighted least squares regression
* QA changes including new format for input object "known" and option to resample known data by site rather than by sample; for bySite option returned results are the average of site-level average statistics; argument order changed for consistancy with calRaster

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
