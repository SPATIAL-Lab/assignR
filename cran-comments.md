## Test environments
* local Windows 10 x64, R 3.6.1
* Ubuntu 16.04.6 (on Travis CI), R 3.6.1
* Windows x64 (on win-builder), R 3.7.0

## R CMD check results
No ERRORs, WARNINGs or NOTEs 

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission changes
* added on.exit() calls to ensure reversion of values of options(scipen) and par(mfrow) set in jointP() and plot.QA(), resectively
* removed dontrun{} from examples in QA.Rd and plot.QA.Rd and stripped down examples to ensure run time < 5s
* added ISBN in place of URL for citation listed in DESCRIPTION
* found and patched backwards compatability issue in QA.R (args to set.seed())