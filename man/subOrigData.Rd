\name{subOrigData}

\alias{subOrigData}

\title{
Extract known origin sample data
}

\description{
This function subsets the known-origin isotope dataset included in this package and conducts optional transformations to convert isotope measurements to a common calibration scale.
}

\usage{
subOrigData(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
  age_code = NULL, mask = NULL, std_scale = "VSMOW_H", niter = 1000)
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
  \item{reference}{
character string or string vector. Bibliographic reference for data to be extracted.
}
  \item{age_code}{
character string or string vector. Animal age code for data to be extracted. 
}
  \item{mask}{
SpatialPolygons or SpatialPolygonsDataFrame. Constrains the geographic area from which data are extracted. If not provided, global.
}
  \item{std_scale}{
character string. Text identifier for standard calibration scale to which all isotope values will be transformed. See \code{\link[assignR]{hsds}} and \code{\link[assignR]{osds}}. Use \dQuote{NULL} to obtain untransformed values on the originally reported standard scale.}
  \item{niter}{
integer. Number of random samples used to propagate uncertainty in standard chain transformations.
}
}

\value{
Returns an object of class \dQuote{subOrigData}, formatted for use in \code{\link[assignR]{calRaster}} or \code{\link[assignR]{QA}} functions.
  \item{data}{
SpatialPointsDataFrame including one feature for each selected sample. data.frame fields are described in \code{\link[assignR]{knownOrig_samples}}. }
  \item{sources}{
data.frame. Information for all data sources for the selected samples. Fields are described in \code{\link[assignR]{knownOrig_sources}} 
  }
  \item{chains}{
list. Each item is a character string vector containing the chain of standard scales used in the transformation for a set of selected samples. See \code{\link[assignR]{hsds}} and \code{\link[assignR]{osds}}.
  }
  \item{marker}{
character string. The isotopic marker specified in the call to \code{subOrigData}
  }
}

\references{
Magozzi, S. et al. (in prep) Standard-chain transformation to improve the comparability of tissue hydrogen and oxygen isotope data.
}

\examples{
## WITHOUT mask
# extract d2H data for Jackdaw, Partridge and Willow Grouse, transformed to the VSMOW/SLAP H standard scale
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
summary(d1)

# extract d2H data for insects and passerine birds without transformation
d2 = subOrigData(group = c("Insect","Passerine"), std_scale = NULL)
summary(d2)

# extract d18O data for all humans, transformed to the VSMOW/SLAP O standard scale
d3 = subOrigData(marker = "d18O", group = c("Modern human", "Indigenous human"), std_scale = "VSMOW_O")
summary(d3)

# extract d2H data for humans using taxon, transformed to the VSMOW/SLAP H standard scale
d4 = subOrigData(marker = "d2H", taxon = "Homo sapiens", std_scale = "VSMOW_H")
summary(d4)

## WITH mask
# error - no samples found
\dontrun{d5 = subOrigData(taxon = "Turdus philomelos", mask = naMap)}
# this works OK
d6 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", 
  "Turdus migratorius"), mask = naMap)
}
