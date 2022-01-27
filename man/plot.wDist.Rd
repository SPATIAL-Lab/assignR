\name{plot.wDist}

\alias{plot.wDist}

\title{
Plot weighted distance and bearing distributions
}

\description{
Plot the output from \code{\link{wDist}}, including weighted kernel density distributions for distance and bearing of travel.
}

\usage{
\method{plot}{wDist}(x, ..., bin = 20, pty = "both", index = c(1:5))
}

\arguments{
  \item{x}{
    A wDist object
  }
  \item{...}{
    Other arguments to be passed to plot
  }
  \item{bin}{
    numeric. Bin width used to generate rose plot of travel bearings, in degrees. Must be a factor of 360.
  }
  \item{pty}{
    character. Type of plot to produce. Must be one of \dQuote{dist}, \dQuote{bear}, or \dQuote{both}.
  }
  \item{index}{
    numeric. Which items in x to plot? Numeric vector of up to 5 integers. Values in excess of 5 or exceeding the length of x will be ignored.
  }
}

\details{
For the default \code{pty}, two plot panels will be printed to the active graphical device showing the distance and bearing distributions for (up to) the first five samples in \code{wd}. If more than five items exist in \code{wd}, those beyond the fifth will be ignored and a message returned.
}

\seealso{
\code{\link{wDist}}
}

\examples{
library(sp)

# load North America boundary and global isoscape
data("naMap")
data("d2h_lrNA")

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

# plot distributions
plot(wd)

# plot bearing distriubtion for sample B with a finer bin size
plot(wd, bin = 5, pty = "bear", index = 2)
}
