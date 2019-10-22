## Test environments
* local Windows 10 x64, R 3.6.1
* Ubuntu 16.04.6 (on Travis CI), R 3.6.1
* Windows x64 (on win-builder), R 3.7.0

## R CMD check results
No ERRORs or WARNINGs 

There was 1 NOTE (win-builder):

Maintainer: 'Gabe Bowen <gabe.bowen@utah.edu>'
New submission

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission changes
* 'Tools for' removed from title, now 'Infer Geographic Origin from Isotopic Data'
* All print()/cat() removed from functions and replaced with message() or warning(): subOrigData.R, oddsRatio.R, and calRaster.R
* Writing to files has been removed from all functions unless user specifies a directory in the function call: calRaster.R, pdRaster.R, qtlRaster.R, plot.QA.R; verified examples and vignettes do not write to disk