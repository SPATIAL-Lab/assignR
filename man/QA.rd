\name{QA}

\alias{QA}

\title{
Quality assessment of geographic assignments
}

\description{
How well does a given isoscape and/or known origin data set constrain the geographic origin of samples? Uses iterative re-sampling of known origin data to evaluate sample assignments and reports a suite of quality metrics.
}

\usage{
QA(known, isoscape, bySite = TRUE, valiStation = 1, valiTime = 50, 
  recal = TRUE, by = 2, prior = NULL, mask = NULL, setSeed = TRUE, 
  name = NULL)
}

\arguments{
  \item{known}{
subOrigData, list of subOrigData, or SpatialPointsDataFrame. Known-origin tissue isotope data from the \code{subOrigData} function or provided by user. User-provided data must be formatted as subOrigData objects (see \code{\link{subOrigData}}) or a SpatialPointsDataFrame (see Details).
}
  \item{isoscape}{SpatRaster with two layers or \code{\link{isoStack}} object. For user-generated raster objects, the first layer must be the isoscape (mean prediction) and the second the isoscape prediction uncertainty (1 standard deviation).
}
  \item{bySite}{logical. Resample known by site (TRUE) or by sample (FALSE)?}
  \item{valiStation}{numeric. How many sites or samples from known are withheld for validation? Must be two or more smaller than the length of \code{known}.
}
  \item{valiTime}{numeric. How many times do you want to randomly draw validation samples and run the validation? Must be an integer greater than one.
}
  \item{recal}{logical. Recalibrate the isoscape(s) using the known-origin data? If FALSE, \code{isoscape} should be a calibrated product appropriate to the samples, and a single iteration is run for each sample in \code{known}; parameters \code{bySite}, \code{valiStation}, and \code{valiTime} are ignored.}
  \item{by}{integer. Threshold increment to use in evaluating assignment performance. Must be between 1 and 25.}
  \item{prior}{raster. Optional raster layer with prior probabilities, which has the same projection, resolution and extent as \code{isoscape}.
}
  \item{mask}{SpatialPolygonsDataFrame. Constrains the area of the output rasters. If this is not provided, the entire area of \code{isoscape} is returned.
}
  \item{setSeed}{logical. Do you want to \code{set.seed()} when you randomly draw validation stations? \dQuote{TRUE} gives the same sequence of random draws each time the function is called.
}
  \item{name}{character. Useful for identifying the QA output in subsequent plotting.
  }
}

\details{
If \code{known} is a user-provided SpatialPointsDataFrame, the first field in \code{@data} must include the measured value for the first (or only) isotope marker and the second the one standard deviation uncertainty on that value. Subsequent fields must include the same information for all other isotope markers included in the analysis, and these markers must appear in the same order as in \code{isoscape}. A user-provided SpatialPointsDataFrame must include a field named \dQuote{Site_ID} containing unique values for each sampling site to support the \dQuote{bySite} option, otherwise use \code{bySite = FALSE}.
}

\value{
Returns an object of class \dQuote{QA}.
\item{val_stations}{numeric. An X*Y data.frame of validation station IDs for all valiTime. X = \code{valiTime} and Y = \code{valiStation}.
}
\item{pd_val}{numeric. An X*Y data.frame containing the posterior probability density for the validation stations. X = \code{valiTime} and Y = \code{valiStation}.
}
  \item{prption_byArea}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 area quantile (with increment of \code{by / 100}; Y = \code{ceiling(100 / by) + 1}). X = \code{valiTime}.
}
  \item{prption_byProb}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 probability quantile (with increment of \code{by / 100}; Y = \code{ceiling(100 / by) + 1}). X = \code{valiTime}.
}
  \item{precision}{list. The length of the list is \code{valiTime}. Each element is an X*Y matrix showing the proportional area of the total assignment surface covered by the assignment region at a given probability quantile from 0.00 to 1.00 (with increment of \code{by / 100}; X = \code{ceiling(100 / by) + 1}) for each validation sample (Y = \code{valiStation}).}
  \item{random_prob_density}{Random probability of assignment to any given grid cell on the assignment surface(i.e. 1 divided by the total number of grid cells).
}
  \item{name}{character. Name assigned to the QA object.}
  \item{by}{integer. Value of by used.}
}

\note{
See Ma et al. (2020) for methodological details.
}

\references{
Ma, C. et al. (2020) assignR : An R package for isotope-based geographic assignment. \emph{Methods in Ecology and Evolution} \strong{11} 996--1001. \doi{10.1111/2041-210X.13426}.

Vander Zanden, H. B. et al. (2014) Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. \emph{Methods in Ecology and Evolution} \strong{5} 891--900. \doi{10.1111/2041-210X.12229}
}

\seealso{
\code{\link[assignR]{plot.QA}}
}

\examples{
library(terra)

# load North America boundary and global isoscape
library(terra)
data("naMap")
d2h_lrNA = rast(system.file("data/d2h_lrNA.tif", package = "assignR"))

# extract some known-origin data
d1 = subOrigData(taxon = "Buteo lagopus")

# run quality assessment based on precipitation hydrogen isotopes and 
# known-origin samples; small values of valiStation and valiTime 
# are used in example to reduce run time

# first with one example
# gives warning because a small number of samples are available
qa1 = QA(known = d1, isoscape = d2h_lrNA, valiTime = 2, by = 10, 
          mask = naMap, name = "Buteo")
                    
# plot the qa result
plot(qa1)

# now compare with a second data set
\donttest{d2 = subOrigData(taxon = "Charadrius montanus")
qa2 = QA(known = d2, isoscape = d2h_lrNA, valiTime = 2, by = 10, 
          mask = naMap, name = "Charadrius")
plot(qa1, qa2)}
}
