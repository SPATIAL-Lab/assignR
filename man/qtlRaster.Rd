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
  RasterLayer, RasterStack, or RasterBrick. Probability density maps for individual samples, e.g., as output by \code{pdRaster}.
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

\examples{
# load North America boundary and global isoscape
data("naMap")
data("d2h_world")

# load hydrogen isotope data for human hair in North America
d = subOrigData(taxon = "Homo sapiens", 
  reference = "Ehleringer et al. 2008", mask = naMap)

# rescale from environmental isoscape to tissue isoscape
r = calRaster(known = d, isoscape = d2h_world, mask = naMap)

# four unknown-origin examples
id = c("A", "B", "C", "D")
d2H = c(-110, -90, -105, -102)
un = data.frame(id, d2H)

# posterior probability surface
asn = pdRaster(r,unknown=un,mask=naMap)

# assign to most probable 10 percent of area
\dontrun{qtlRaster(asn, threshold = 0.1, savePDF = TRUE)}

# assign to most probable 10 percent of proabability distribution
qtlRaster(asn, threshold = 0.1, thresholdType = "prob")
}
