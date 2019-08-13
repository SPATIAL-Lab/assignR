\name{plot.QA}
\alias{plot.QA}
\title{
PLot the quality assessment of the geographic assignment.
}
\description{
Plot the output from QA results. Accuracy, precision and probability density of validation stations.
}
\usage{
plot.QA(obj, savePDF = FALSE)
}
\arguments{
  \item{obj}{
    Either a single QA object or a list of QA objects
  }
  \item{savePDF}{
    logical. Save the plot as PDF (TRUE) or not (FALSE [default]).
  }
}
\references{
Ma et al., 2019

Vander Zanden, H.B., Wunder, M.B., Hobson, K.A., Van Wilgenburg, S.L., Wassenaar, L.I., Welker, J.M. and Bowen, G.J., 2014. Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. Methods in Ecology and Evolution, 5(9), pp.891-900.
}

\seealso{
\code{\link[ggplot2]{ggplot}}
\code{\link[assignR]{QA}}
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
bird_d2h <- bird[1:20,c("Longitude", "Latitude", "d2H")]
coordinates(bird_d2h) <- c(1,2)
proj4string(bird_d2h) <- proj4string(d2h_world)

# run quality assessment based hydrogen isotope from precipitation and known-origin bird
d2h_QA <- QA(isoscape = r, known = bird_d2h, valiStation = 2, 
                    valiTime = 5, setSeed = T)
                    
# plot the QA result
plot(d2h_QA)
}
