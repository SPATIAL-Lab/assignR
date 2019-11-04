## Test environments
* local Windows 10 x64, R 3.6.1
* Ubuntu 16.04.6 (on Travis CI), R 3.6.1
* Windows x64 (on win-builder), R 3.7.0

## R CMD check results
No ERRORs, WARNINGs or NOTEs 

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission changes
* Reduced run-time for vignette by reducing raster resolution and resampling iterations. win-builder (release version) now showing vignette rebuild time of 215s and total build time of 537s.