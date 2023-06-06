\name{wDist}

\alias{wDist}

\title{
Probability weighted distances and bearings
}

\description{
Calculate the distance and bearing of migration for one or more samples, weighted by probabilities from a \code{pdRaster} analysis. 
}

\usage{
wDist(pdR, sites, maxpts = 1e5)
}

\arguments{
  \item{pdR}{
  SpatRaster of n probability density maps, e.g., as produced by \code{pdRaster}.
  }
  \item{sites}{
  SpatVector object containing the collection locations for the n samples represented in \code{pdR}.
  }
  \item{maxpts}{
  numeric. Maximum number of grid cells at which to calculate bearing and distance.
  }
}

\details{
\code{pdR} and \code{sites} must be of equal length and corresponding order. Names in the output object are taken from the names of the layers in \code{pdR}.

Distances and bearings are calculated on the WGS84 geoid using functions from the \pkg{terra} and \pkg{geosphere} package. These calculations can take a long time for large rasters. If \code{maxpts} is less than the number of grid cells in each \code{pdR} layer, calculations are carried out for \code{maxpts} randomly selected cells.

Bearing values correspond to the initial bearing from source to collection location, and are reported on a scale of -180 to +180 degrees. The statistical metrics are rectified so that values for distributions spanning due south are reported correctly. Both weighted bearing and distance distributions are often multimodal, and it is recommended to review the distribution densities to assess the representativeness of the statistics (e.g., using \code{\link{plot.wDist}}). 
}

\value{
Returns an object of class \dQuote{wDist}, a list of length n. Each item contains three named objects:
  \item{stats}{
named number. Statistics for the distance (dist, meters) and bearing (bear, degrees) between source and collection locations, including the weighted mean (wMean) and quantile (wXX) values.}
  \item{d.dens}{
density. Weighted kernel density for the distance between source and collection locations (meters). See \code{\link[stats]{density}}. 
  }
  \item{b.dens}{
density. Weighted kernel density for the bearing between source and collection locations (degrees). See \code{\link[stats]{density}}.
  }
}

\examples{
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

# random collection locations
sites = d$data[sample(seq(length(d$data)), 4),]

# generate a wDist object
wd = wDist(pp, sites)

# structure of the wDist object
str(wd, 2)
}

