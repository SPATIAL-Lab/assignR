\name{subOrigData}
\alias{subOrigData}
\title{
subOrigData function
}
\description{
This function subsets the known-origin isotope dataset included in this package
}
\usage{
subOrigData(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
  age_code = NULL, mask = NULL)
}
\arguments{
  \item{marker}{
String. Column name for isotopic data to be extracted. Possible values include "d2H" and "d18O".
}
  \item{taxon}{
String or a string vector. Species names for data to be extracted. Details can be seen in \code{\link[assignR]{knownOrig}}
}
  \item{group}{
String or a string vector. Organism groups for data to be extracted. Details can be seen in \code{\link[assignR]{knownOrig}}
}
  \item{reference}{
String or a string vector. Bibliographic reference for data to be extracted. Details can be seen in \code{\link[assignR]{knownOrig}}
}
  \item{age_code}{
String or a string vector. Animal age code for data to be extracted. Details can be seen in \code{\link[assignR]{knownOrig}}
}
  \item{mask}{
SpatialPolygons or SpatialPolygonsDataFrame. Constrains the geographic area from which data are extracted. If not provided, default of whole world is used.
}
}
\value{
Data selected will be returned as SpatialPointsDataFrame with a single feature corresponding to the selected marker. Output is formatted for use in calRaster or QA functions.
}


\examples{
## WITHOUT mask
# extract d2H data for Jackdaw, Partridge and Willow Grouse
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
summary(d1)

# extract d2H data for [1]insects and [2]birds
d2 = subOrigData(group = c("insect","bird"))
summary(d2)

# extract d18O data for humans
d3 = subOrigData(marker = "d18O", group = "human")
summary(d3)

# extract d2H data for humans using taxon
d4 = subOrigData(marker = "d2H", taxon = "Homo sapiens")
summary(d4)

## WITH mask
d5 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"), 
  mask = naMap)
d6 = subOrigData(group = "human", mask = naMap)
}
