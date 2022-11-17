\name{plot.isoStack}

\alias{plot.isoStack}

\title{
Plot method for stacked isoscapes
}

\description{
Plot the output from \code{\link{isoStack}}.
}

\usage{
\method{plot}{isoStack}(x, ...)
}

\arguments{
  \item{x}{
    An isoStack object
  }
  \item{...}{
   Other arguments to be passed to plot
  }
}

\seealso{
\code{\link{isoStack}}
}

\examples{
library(terra)

# load two isoscapes
d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))
sr_MI = rast(system.file("data/sr_MI.tif", package = "assignR"))

#stack them
h_s = isoStack(d2h_lrNA, sr_MI)

#plot isoStack
plot(h_s)
}
