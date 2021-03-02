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
data("d18o_world")
data("d2h_lrNA")

#stack them
o_h = isoStack(d18o_world, d2h_lrNA)

#plot isoStack
plot(o_h)
}
