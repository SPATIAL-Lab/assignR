\name{c.wDist}

\alias{c.wDist}

\title{
Combine method for wDist objects
}

\description{
Combine statistics from one or more \code{\link{wDist}} objects in a single data frame.
}

\usage{
\method{c}{wDist}(...)
}

\arguments{
  \item{...}{
    One or more wDist objects
  }
}

\value{
data.frame containing sample IDs, distance, and bearing statistics for each sample in \code{...}
}

\seealso{
\code{\link{wDist}}
}

\examples{
# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, niter = 100)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id,d2H)

# posterior probabilities
pp = pdRaster(r, unknown = un, mask = naMap)

# random collection locations
sites = d$data[sample(seq(length(d$data)), 4),]

# generate a wDist object
wd = wDist(pp, sites)

# combine stats and print
c(wd)
}
