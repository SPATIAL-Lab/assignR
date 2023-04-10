v2.2.3 fixes a new test error on Mac arm64.

# Submission comment
assignR was archived despite the triggering error having been fixed in v2.2.2. Apparently due to the new test error that is fixed here (and emerged after testing/submission/acceptance of v2.2.2).

# Test environments
* Mac OS 11.5.2 (on macbuilder); R 4.2.1
* local Windows 11 x64; R 4.2.3
* Ubuntu 22.04.2 (on GitHub Actions); R 4.2.3
* Mac OS 12.6.3 (on GitHub Actions); R 4.2.3
* Windows Server 2022 (on GitHub Actions); R 3.6
* Windows Server 2022 (on GitHub Actions); r84036

# R CMD check results
No ERRORs, WARNINGs, or NOTEs

# Downstream dependencies
There are currently no downstream dependencies for this package.
