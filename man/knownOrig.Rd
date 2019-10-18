\name{knownOrig}

\alias{knownOrig}

\docType{data}

\title{
Isotope values for tissues of known origin
}

\description{
This dataset consists of hydrogen and oxygen isotope values of human hair, insect chitin and bird feathers with location and sample information.
}

\usage{data("knownOrig")
}

\format{
SpatialPointsDataFrame with 8 variables.
\describe{
  \item{[, 1]}{d2H: Hydrogen isotope ratio (permil, VSMOW)}
  \item{[, 2]}{d18O: Oxygen isotope ratio (permil, VSMOW)}
  \item{[, 3]}{Taxon: Genus and species name}
  \item{[, 4]}{Group: Taxonomic group (informal)}
  \item{[, 5]}{Source_quality:Code indicating level of certainty in geographic origin}
  \item{[, 6]}{Age_code: Code for age of individual}
  \item{[, 7]}{Reference: Bibliographic citation of source}
  \item{[, 8]}{ID: Unique ID}
}
}

\references{
Bowen, G. J. et al. (2009) Dietary and physiological controls on the hydrogen and oxygen isotope ratios of hair from mid-20th century indigenous populations. \emph{American Journal of Physical Anthropology}, \strong{139}, 494--504.

Clark, R. G., Hobson, K. A. and Wassenaar, L. I. (2006) Geographic variation in the isotopic (dD, d13C, d15N, d34S) composition of feathers and claws from lesser scaup and northern pintail: implications for studies of migratory connectivity. \emph{Canadian Journal of Zoology} \strong{84} 1395--1401.

Ehleringer, J. R. et al. (2008) Hydrogen and oxygen isotope ratios in human hair are related to geography. \emph{Proceedings of the National Academy of Sciences} \strong{105} 2788--2793.

Hobson, K. A., Bowen, G. J., Wassenaar, L. I., Ferrand, Y. and Lormee, H. (2004) Using stable hydrogen and oxygen isotope measurements of feathers to infer geographical origins of migrating European birds. \emph{Oecologia} \strong{141} 477--488.

Hobson, K. A. and Koehler, G. (2015) On the use of stable oxygen isotope (d18O) measurements for tracking avian movements in North America. \emph{Ecology and Evolution} \strong{5} 799--806.

Hobson, K. A., Van Wilgenburg, S. L., Wassenaar, L. I. and Larson, K. (2012) Linking hydrogen (d2H) isotopes in feathers and precipitation: sources of variance and consequences for assignment to isoscapes. \emph{PLoS ONE} \strong{7} e35137.

Hobson, K. A. and Wassenaar, L. I. (1997) Linking breeding and wintering grounds of neotropical migrant songbirds using stable hydrogen isotopic analysis of feathers. \emph{Oecologia} \strong{109} 142--148.

Hobson, K. A., Wassenaar, L. I. and Taylor, O. R. (1999) Stable isotopes (dD and d13C) are geographic indicators of natal origins of monarch butterflies in easthern North America. \emph{Oecologia} \strong{120} 397--404.

Lott, C. A. and Smith, J. P. (2006) A geographic-information-system approach to estimating the origin of migratory raptors in North America using stable hydrogen isotope ratios in feathers. \emph{The Auk} \strong{123} 822--835.

Thompson, A. H. et al. (2010) Stable isotope analysis of modern human hair collected from Asia (China, India, Mongolia and Pakistan). \emph{American Journal of Physical Anthropology} \strong{141} 440--451.

Wunder, M.B., 2007. Geographic structure and dynamics in Mountain Plover. Doctoral dissertation, Colorado State University.
}
\seealso{
\code{\link[assignR]{subOrigData}}
}
\examples{
data("knownOrig")
data("wrld_simpl")
plot(wrld_simpl)
plot(knownOrig, add = TRUE, col = "red", cex = 0.5)
}
\keyword{datasets}
