## Test environments
* local Windows 10 x64, R 3.6.3
* Ubuntu 16.04.6 (on Travis CI), R 4.0.0
* Windows x64, i386 (on win-builder), R 4.0.0
* Windows x64, i386 (on win-builder), r78617

## R CMD check results
No ERRORs, or WARNINGs

1 NOTE (win-builder current release only):

Examples with CPU (user + system) or elapsed time > 10s
   user system elapsed
QA 9.31   0.61   10.25

This example clocks at right around 10s on win-builder, usually a little less. We've worked to reduce run-time w/ each new package version and keep this at or under 10 seconds; the current version demonstrates the minimum functionality. 

## Downstream dependencies
There are currently no downstream dependencies for this package
