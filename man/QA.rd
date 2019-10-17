\name{QA}
\alias{QA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the efficacy of geographic assignment
}
\description{
What is the power of a certain isoscape used for geographic assignment? Using the known origin data and the isoscape as input to test it. You will get the population accuracy, precision and probability density (see returned ).
}
\usage{
QA(isoscape, known, valiStation = floor(length(known)*0.1), valiTime = 50, 
  mask = NULL, setSeed = TRUE, name = NULL)
}
\arguments{
  \item{isoscape}{raster. Environmental isoscape. Two layers: the first one is mean and the second one is standard deviation.}
  \item{known}{SpatialPointsDataFrame. known-origin data that should contain one feature: tissue isotope value that should be the same isotope as the environmental isoscape. This input must have a coordinate reference system as that in isoscape.
}
  \item{valiStation}{numeric. How many stations of the known origin with tissue isotope are used for validation. This must be smaller than the total number of known origin.
}
  \item{valiTime}{numeric. How many times do you want to randomly draw validation stations and run the validation. Must be an integer greater than one. 
}
  \item{mask}{SpatialPolygonsDataFrame. Constrains the area of the output rasters. If this is not provided, the entire area of isoscape is returned.
}
  \item{setSeed}{logical. Do you want to set.seed() when you randomly draw validation stations? Yes gives the same sequence of random draws each time the function is called.
}
  \item{name}{character. Useful for identifying the QA output in subsequent plotting.
  }
}

\value{
\item{val_stations}{numeric. An X*Y data.frame of validation station IDs for all valiTime. X = valiTime and Y contains the validation station row IDs.}
\item{pd_val}{numeric. An X*Y data.frame containing the posterior probability density for the validation stations. X = valiTime and Y = valiStation.
}
  \item{prption_byArea}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 area quantile (with increment of 0.01; Y = 101). X = valiTime.
}
  \item{prption_byProb}{numeric. An X*Y data.frame showing the proportion of validation individuals for which the known origin is contained within the top 0.00 to 1.00 probability quantile (with increment of 0.01; Y = 101). X = valiTime.
}
  \item{precision}{list. The length of the list is valiTime. Each element is an X*Y matrix showing the proportional area of the total assignment surface covered by the assignment region at a given probability quantile from 0.00 to 1.00 *with increment of 0.01; X = 101) for each validation sample (Y = valiStation).}
  \item{random_prob_density}{Random probability of assignment to any given gridcell on the assignment surface(i.e. 1 divided by the total number of grid cells).
}
  \item{name}{character. Name assigned to the QA object.}
}

\note{
Please see Ma et al., 2019 for details of these values returned and the methodology.
}

\references{
Vander Zanden, H.B., Wunder, M.B., Hobson, K.A., Van Wilgenburg, S.L., Wassenaar, L.I., Welker, J.M. and Bowen, G.J., 2014. Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. Methods in Ecology and Evolution, 5(9), pp.891-900.
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
d1 = subOrigData(taxon = "Lanius ludovicianus")
d2 = subOrigData(taxon = "Buteo lagopus")

\dontrun{# run quality assessment based hydrogen isotope from precipitation and known-origin bird
qa1 = QA(isoscape = d2h_world, known = d1, valiStation = 2, 
                    valiTime = 5, mask = naMap, name = "Lanius")
                    
qa2 = QA(isoscape = d2h_world, known = d2, valiStation = 2, 
                    valiTime = 5, mask = naMap, name = "Buteo")

# plot the QA result
plot(qa1, qa2)}
}
