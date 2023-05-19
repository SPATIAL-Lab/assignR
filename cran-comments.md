v2.2.3 fixes a new test error on Mac arm64.

# Submission comment
v2.2.2 was archived. The problem triggering the archive countdown was already fixed in v2.2.2. A new platform-specific test error emerged after testing/submission/acceptance of v2.2.2 and has been fixed here and verified on macbuilder.

# Test environments
* Mac OS 11.5.2 (on macbuilder); R 4.2.1
* local Windows 11 x64; R 4.2.3
* Ubuntu 22.04.2 (on GitHub Actions); R 4.2.3
* Mac OS 12.6.4 (on GitHub Actions); R 4.2.3
* Windows Server 2022 (on GitHub Actions); R 3.6
* Windows Server 2022 (on GitHub Actions); r84210

# R CMD check results
No ERRORs, WARNINGs, or NOTEs

# Downstream dependencies
There are currently no downstream dependencies for this package.
