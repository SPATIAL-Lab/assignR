\name{subOrigData}

\alias{subOrigData}

\title{
Extract known origin sample data
}

\description{
This function subsets the known-origin isotope dataset included in this package.
}

\usage{
subOrigData(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
  age_code = NULL, mask = NULL)
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
SpatialPolygons or SpatialPolygonsDataFrame. Constrains the geographic area from which data are extracted. If not provided, default of global is used.
}
}

\value{
SpatialPointsDataFrame including the location and a single feature corresponding to the marker for each selected sample. Output is formatted for use in calRaster or QA functions.
}

\examples{
## WITHOUT mask
# extract d2H data for Jackdaw, Partridge and Willow Grouse
d1 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", "Turdus migratorius"))
summary(d1)

# extract d2H data for insects and birds
d2 = subOrigData(group = c("insect","bird"))
summary(d2)

# extract d18O data for humans
d3 = subOrigData(marker = "d18O", group = "human")
summary(d3)

# extract d2H data for humans using taxon
d4 = subOrigData(marker = "d2H", taxon = "Homo sapiens")
summary(d4)

## WITH mask
# error - no samples found
\dontrun{d5 = subOrigData(taxon = "Turdus philomelos", mask = naMap)}
# these work OK
d6 = subOrigData(taxon = c("Danaus plexippus", "Setophaga ruticilla", 
  "Turdus migratorius"), mask = naMap)
d7 = subOrigData(group = "human", mask = naMap)
}
