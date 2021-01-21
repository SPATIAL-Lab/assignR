\name{isoStack}

\alias{isoStack}

\title{
Stack isoscapes
}

\description{
Combine multiple isoscapes into a single data object, including optional reconsiliation of raster properties.
}

\usage{
isoStack(..., clean = TRUE)
}

\arguments{
  \item{...}{
  Two or more RasterStack or RasterBrick isoscapes, each with two layers, e.g., as produced by \code{\link{calRaster}}. The first layer in each should contain the tissue-specific isoscape and the second the isoscape prediction uncertainty (1 standard deviation).
  }
  \item{clean}{
  logical. Reconcile differences in raster properties within \code{...}?
  }
}

\details{
If \code{clean} = TRUE all rasters are projected to the projection of the first object in \code{...} and then resampled to the highest spatial resoluction and smallest common spatial extent within \code{...}. Finally, cells containing NA in any layer within \code{...} are masked across all layers.

If \code{clean} = FALSE any differences in raster properties between isoscapes will produce an error. 
}

\value{
Returns an object of class \dQuote{isoStack}, a list containing the isoscapes objects in \code{...} after any cleaning. 
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