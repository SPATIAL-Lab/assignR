\name{qtlRaster}

\alias{qtlRaster}

\title{
Sample assignment using thresholds
}

\description{
Selects the grid cells of probability density rasters with the highest probability and returns rasters with these cell values set to 1. Cells are selected based on the user-specified quantile threshold so that the most-probable cells representing a given fraction of the assignment area or posterior probability are returned.
}

\usage{
qtlRaster(pdR, threshold, thresholdType = "area", genplot = TRUE, outDir = NULL)
}

\arguments{
  \item{pdR}{
  SpatRaster. Probability density maps for individual samples, e.g., as output by \code{pdRaster}.
  }
  \item{threshold}{
  numeric from 0 to 1. Quantile to be selected. 
  }
  \item{thresholdType}{
  character. Either \dQuote{area} (default) or \dQuote{prob}. If \dQuote{area}, the most probable cells constituting \code{threshold} percent of the assignment area are selected. If \dQuote{prob}, the most probable cells constituting \code{threshold} percent of the posterior probability are selected.
  }
  \item{genplot}{
  logical.Plot results in R.
  }
  \item{outDir}{
  character string. Directory path to which output will be saved. If NULL no files are written.
  }
}

\value{
SpatRaster including a binary assignment surface for each individual in \code{pdR}. If \code{outDir} is not NULL, writes individual rasters in GeoTIFF format and a single PDF file with images for each raster to the designated directory.
}

\examples{
library(terra)

# load hydrogen isotope data for human hair in North America
d = subOrigData(group = "Modern human", mask = naMap, genplot = FALSE)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(d, d2h_lrNA, naMap, genplot = FALSE)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id, d2H)

# posterior probabilities
pp = pdRaster(r, un, mask = naMap, genplot = FALSE)

# assign to most probable 10 percent of area
\dontrun{qtlRaster(pp, threshold = 0.1)}

# assign to most probable 10 percent of proabability distribution
qtlRaster(pp, threshold = 0.1, thresholdType = "prob")
}
