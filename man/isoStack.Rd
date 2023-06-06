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
  Two or more SpatRaster isoscapes, each with two layers, or \code{rescale} objects as produced by \code{\link{calRaster}}. For objects other than \code{rescale}, the first layer in each should contain the tissue-specific isoscape and the second the isoscape prediction uncertainty (1 standard deviation).
  }
  \item{clean}{
  logical. Reconcile differences in raster properties within \code{...}?
  }
}

\details{
If \code{clean} = TRUE all raster layers are projected to the projection of the first object in \code{...} and then resampled to the highest spatial resolution and smallest common spatial extent within \code{...}. Finally, cells containing NA in any layer within \code{...} are masked across all layers.

If \code{clean} = FALSE any differences in raster properties between isoscapes will produce an error. 
}

\value{
Returns an object of class \dQuote{isoStack}, a list containing the isoscapes objects in \code{...} after any cleaning. 
}

\examples{
#stack H and Sr isoscapes
h_s = isoStack(d2h_lrNA, sr_MI)
}
