\name{QA}

\alias{QA}

\title{
Quality assessment of geographic assignments
}

\description{
How well does a given isoscape and known origin data set constrain the geographic origin of samples? Uses iterative re-sampling of the known origin data set to evaluate sample assignments and reports a suite of quality metrics.
}

\usage{
QA(isoscape, known, valiStation = ceiling(length(known)*0.1), valiTime = 50, 
  mask = NULL, setSeed = TRUE, name = NULL)
}

\arguments{
  \item{isoscape}{RasterStack or RasterBrick with two layers, e.g., as created by \code{calRaster}. The first layer is the tissue-specific isoscape and the second the isoscape prediction uncertainty (1 standard deviation).
}
  \item{known}{SpatialPointsDataFrame. Known-origin data that should contain only one feature: tissue isotope value. Its length must be larger or equal to 3. Known-origin data can be queried using \code{knownOrig}.
}
  \item{valiStation}{numeric. How many samples from known are withheld for validation? Must be two or more smaller than the length of \code{known}.
}
  \item{valiTime}{numeric. How many times do you want to randomly draw validation stations and run the validation? Must be an integer equal to or greater than one. 
}
  \item{mask}{SpatialPolygonsDataFrame. Constrains the area of the output rasters. If this is not provided, the entire area of \code{isoscape} is returned.
}
  \item{setSeed}{logical. Do you want to \code{set.seed()} when you randomly draw validation stations? Yes gives the same sequence of random draws each time the function is called.
}
  \item{name}{character. Useful for identifying the QA output in subsequent plotting.
  }
}

\value{
\item{val_stations}{numeric. An X*Y data.frame of validation station IDs for all valiTime. X = \code{valiTime} and Y = \code{valiStation}.
}
\item{pd_val}{numeric. An X*Y data.frame containing the posterior probability density for the validation stations. X = \code{valiTime} and Y = \code{valiStation}.
}
  \item{prption_byArea}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 area quantile (with increment of 0.01; Y = 101). X = \code{valiTime}.
}
  \item{prption_byProb}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 probability quantile (with increment of 0.01; Y = 101). X = \code{valiTime}.
}
  \item{precision}{list. The length of the list is \code{valiTime}. Each element is an X*Y matrix showing the proportional area of the total assignment surface covered by the assignment region at a given probability quantile from 0.00 to 1.00 *with increment of 0.01; X = 101) for each validation sample (Y = \code{valiStation}).}
  \item{random_prob_density}{Random probability of assignment to any given grid cell on the assignment surface(i.e. 1 divided by the total number of grid cells).
}
  \item{name}{character. Name assigned to the QA object.}
}

\note{
Please see Ma et al., 2019 for methodological details.
}

\references{
Ma et al. (in review) Does transpiration matter? Comparing geographic assignment with precipitation- and plant-based isoscapes using IsoMAP and assignR software. \emph{Movement Ecology}.

Vander Zanden, H. B. et. al (2014) Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. \emph{Methods in Ecology and Evolution} \strong{5} 891--900.
}

\seealso{
\code{\link[assignR]{plot.QA}}
}

\examples{
# load data
data("naMap") # North America 
data("d2h_world") # precipitation hydrogen isotope of the world
data("knownOrig") # hydrogen isotopes of known-origin samples

# extract some known-origin data
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")

# run quality assessment based on precipitation hydrogen isotopes and 
# known-origin birds; small values of valiStation and valiTime 
# used in example to reduce run time

# first with one example
qa1 = QA(isoscape = d2h_lrNA, known = d1, valiStation = 1, 
          valiTime = 2, mask = naMap, name = "Charadrius")
                    
# plot the qa result
plot(qa1)

# now compare
\donttest{qa2 = QA(isoscape = d2h_lrNA, known = d2, valiStation = 1, 
          valiTime = 2, mask = naMap, name = "Buteo")
plot(qa1, qa2)}
}
