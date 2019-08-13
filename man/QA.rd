\name{QA}
\alias{QA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the efficacy of geographic assignment using a certain type of isoscape
}
\description{
What is the power of a certain isoscape used for geographic assignment? Using the known origin data and the isoscape as input to test it. You will get the population accuracy, precision and probability density (see returned ).
}
\usage{
QA(isoscape, known, valiStation, valiTime, setSeed = T)
}
\arguments{
  \item{isoscape}{raster. Environmental isoscape. Two layers: the first one is mean and the second one is standard deviation.}
  \item{known}{
SpatialPointsDataFrame. known-origin data that should contain one feature: tissue isotope value that should be the same isotope as the environmental isoscape. This input must have a coordinate reference system as that in isoscape
}
  \item{valiStation}{
  numeric. How many stations of the known origin with tissue isotope are used for validation. This must be smaller than the total number of known origin.
}
  \item{valiTime}{
numeric. How many times do you want to randomly draw validation stations and run the validation.
}
  \item{setSeed}{
Do you want to set.seed() when you run the randomly draw validation stations? If yes and your input data are the same, the output would be exactely the same.
}
}

\value{
\item{val_stations}{numeric. An X*Y data.frame of validation station IDs for all valiTime. X is the valiTime and Y is the validation station IDs.}
\item{pd_bird_val}{
numeric. An X*Y data.frame containing the posterior probability density for the validation stations. X is the simulation numbers = valiTime and Y is the total number of validation stations (valiStation).
}
  \item{prption_byArea}{
  numeric. An X*Y data.frame shows the population-level accuracy that is measured as the proportion of validation individuals in which the known origin is contained within the top probability ranging from 0.01 to 0.99 with increment of 0.01 (99 probabilities which is Y). X is the valiTime. Higher proportion with lower top probability means higher accuracy.
}
  \item{prption_byProb}{
numeric. An X*Y data.frame shows the population-level accuracy that is measured as the proportion of validation individuals in which the known origin is contained within the top percent of area ranging from 0.01 to 0.99 with increment of 0.01 (99 percentage which is Y). X is the valiTime.
}
  \item{precision}{
list. The length of the list is valiTime which reprsents the population-level precision for each validation. It is assessed by the areal proportion of the total surface area covered by the assignment region for each top percent of probability density (threshold). Thresholds from 0.01 to 0.99 with increment of 0.01 are tested. Each elements of the list is a numerical data.frame with dimention of X*Y. X is the 99 threshold. The length of Y is valiStation. Lower proportion with lower threshold in this assessment suggests higher precision.
}
  \item{random_prob_density}{
   Random probability density based on the size of the environmental isoscape (i.e. 1 divided by the total number of grid cells of the isoscape).
}
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
data(naMap) # North America
data(d2h_world) # precipitation hydrogen isotope of the world
data(bird_isotope) # oxygen and hydrogen isotopes of known-origin bird

# crop the world hydrogen data to North America
r <- crop(d2h_world, naMap)
plot(r)

# convert 2 standard deviation from d2h_world to 1 standard deviation
r[[2]] <- r[[2]]/2

# seperate the hydrogen isotope for the known-origin bird
bird_d2h <- bird_isotope[1:20,c("Longitude", "Latitude", "d2H")]
coordinates(bird_d2h) <- c(1,2)
proj4string(bird_d2h) <- proj4string(d2h_world)

# run quality assessment based hydrogen isotope from precipitation and known-origin bird
d2h_QA <- QA(isoscape = r, known = bird_d2h, valiStation = 2,
                    valiTime = 5, setSeed = T)

# plot the QA result
plot(d2h_QA)
}
