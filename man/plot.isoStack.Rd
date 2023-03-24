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
# load two isoscapes
data("d2h_lrNA")
data("sr_MI")

#stack them
h_s = isoStack(list(d2h_lrNA, sr_MI))

#plot isoStack
plot(h_s)
}
