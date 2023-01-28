\name{subOrigData}

\alias{subOrigData}

\title{
Extract known origin sample data
}

\description{
This function subsets the known-origin isotope dataset included in this package and conducts optional transformations to convert isotope measurements to a common reference scale.
}

\usage{
subOrigData(marker = "d2H", taxon = NULL, group = NULL, dataset = NULL, 
  age_code = NULL, mask = NULL, ref_scale = "VSMOW_H", niter = 5000, genplot = TRUE)
}

\arguments{
  \item{marker}{
character string. Column name for isotopic data to be extracted, either \dQuote{d2H} or \dQuote{d18O}.
}
  \item{taxon}{
character string or string vector. Species name(s) for data to be extracted.
}
  \item{group}{
character string or string vector. Taxonomic groups for data to be extracted.
}
  \item{dataset}{
integer or integer vector. Dataset_ID(s) for data to be extracted. See \code{\link[assignR]{knownOrig}} feature \emph{sources}.
}
  \item{age_code}{
character string or string vector. Animal age code for data to be extracted. 
}
  \item{mask}{
SpatialPolygons or SpatialPolygonsDataFrame. Constrains the geographic area from which data are extracted. If not provided, global.
}
  \item{ref_scale}{
character string. Text identifier for reference scale to which all isotope values will be transformed. See \code{\link[assignR]{stds}}. Use \dQuote{NULL} to obtain untransformed values on the originally reported scale.}
  \item{niter}{
integer. Number of random samples used to propagate uncertainty in calibration hierarchy transformations.
}
\item{genplot}{
logical. Plot results in R.
}
}

\value{
Returns an object of class \dQuote{subOrigData}, formatted for use in \code{\link{calRaster}} or \code{\link{QA}} functions.
  \item{data}{
SpatialPointsDataFrame including one feature for each selected sample. data.frame fields are described in \code{\link{knownOrig}} feature \emph{samples}.}
  \item{sources}{
data.frame. Information for all data sources for the selected samples. Fields are described in \code{\link{knownOrig}} feature \emph{sources} 
  }
  \item{chains}{
list. Each item is a character string vector containing the hierarchy of calibrations used in the transformation for a set of selected samples. See \code{\link{stds}}.
  }
  \item{marker}{
character string. The isotopic marker specified in the call to \code{subOrigData}
  }
}

\references{
Magozzi, S. et al. (in press) Calibration chain transformation to improve the comparability of organic hydrogen and oxygen isotope data. \emph{Methods in Ecology and Evolution}
}

\examples{
## WITHOUT mask
# extract d2H data for Jackdaw, Partridge and Willow Grouse, transformed
# to the VSMOW/SLAP H reference scale by default
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", 
  "Turdus migratorius"), niter = 100)
summary(d1)

# extract d2H data for insects and passerine birds without transformation
d2 = subOrigData(group = c("Insect","Passerine"), ref_scale = NULL, genplot = FALSE)
summary(d2)

# extract d18O data for all humans, transformed to the VSMOW/SLAP O reference scale
d3 = subOrigData(marker = "d18O", 
  group = c("Modern human", "Indigenous human"), ref_scale = "VSMOW_O", 
  niter = 100, genplot = FALSE)
summary(d3)

# extract d2H data for humans using taxon, transformed to the VSMOW/SLAP H reference scale
d4 = subOrigData(marker = "d2H", taxon = "Homo sapiens", ref_scale = 
  "VSMOW_H", niter = 100, genplot = FALSE)
summary(d4)

## WITH mask
# error - no samples found
\dontrun{d5 = subOrigData(taxon = "Turdus philomelos", mask = naMap)}
# this works OK
d6 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", 
  "Turdus migratorius"), mask = naMap, genplot = FALSE)
}
