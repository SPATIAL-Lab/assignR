\name{isoStack}

\alias{isoStack}

\title{
Stack isoscapes
}

\description{
Combine multiple isoscapes into a single data object, including optional reconciliation of raster properties.
}

\usage{
isoStack(..., clean = TRUE)
}

\arguments{
  \item{...}{
  Two or more RasterStack or RasterBrick isoscapes, each with two layers, or \code{rescale} objects as produced by \code{\link{calRaster}}. For objects other than \code{rescale}, the first layer in each should contain the tissue-specific isoscape and the second the isoscape prediction uncertainty (1 standard deviation).
  }
  \item{clean}{
  logical. Reconcile differences in raster properties within \code{...}?
  }
}

\details{
If \code{clean} = TRUE all rasters are projected to the projection of the first object in \code{...} and then resampled to the highest spatial resolution and smallest common spatial extent within \code{...}. Finally, cells containing NA in any layer within \code{...} are masked across all layers.

If \code{clean} = FALSE any differences in raster properties between isoscapes will produce an error. 
}

\value{
Returns an object of class \dQuote{isoStack}, a list containing the isoscapes objects in \code{...} after any cleaning. 
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
